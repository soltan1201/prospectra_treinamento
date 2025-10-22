// definir parametros
var vis = {
    mosaico: {
        min: 0.03, max: 0.16,
        bands: ["red","green","blue"]
    },
    vis_pc:{
        min: -0.49, max: -0.24,
    },
    layer_PC: {
        min: -0.5, max: -0.24,
        palette: ["#f09e51", "#f0cb58", "#a8e4a0", "#0e630e"]
    }
}

var param = {
    asset_landsat: "LANDSAT/LC08/C02/T1_RT_TOA",
    asset_id: "LANDSAT/LC08/C02/T1_RT_TOA/LC08_217071_20230713",
    asset_municipiosBA: "users/CartasSol/shapes/municipios_ba",
    bands_names: ["blue", "green", "red", "nir", "swir1", "siwr2"],
    bandas_PCA: ['pc1', 'pc2', 'pc3', 'pc4', 'pc5', 'pc6'],
    lst_bnds : ["B2","B3","B4","B5","B6","B7"],
};  

// funcão 
function calcula_PC(arrayImagem){
    // região de interesse
    var my_geom = ee.Geometry(ee.Image(arrayImagem).geometry());
    // Transforme a imagem em um array
    arrayImagem = arrayImagem.toArray();

    // Calcule a matriz de covariância
    var covariancia = arrayImagem.reduceRegion({
            reducer: ee.Reducer.centeredCovariance(),
            geometry: my_geom,
            scale: 30,
            maxPixels: 1e13
    });

    // Obtenha os autovetores e autovalores
    var covArray = ee.Array(covariancia.get('array'));
    var eigens = covArray.eigen();
    var autovetores = eigens.slice(1, 1);

    // Transforme a imagem original
    var imagemArray2D = arrayImagem.toArray(1);
    var rasterPCA = ee.Image(autovetores).matrixMultiply(imagemArray2D)
                        .arrayProject([0])
                        .arrayFlatten([param.bandas_PCA]);

    return rasterPCA;
}

//  parte principal /////////////////

// Carregue sua imagem
var imagem = ee.Image(param.asset_id).select(param.lst_bnds).rename(param.bands_names);
print("mostrar metadados ", imagem);

//  calculo PCA
var componentesPrincipais = calcula_PC(imagem);
// imgCol.map(calcula_PC)

Map.addLayer(imagem, vis.mosaico, "imagem", true);
// Visualize a primeira componente principal
// Map.addLayer(componentesPrincipais.select('pc1'), vis.vis_pc, 'PC1 simple');
// Map.addLayer(componentesPrincipais.select('pc1'), vis.layer_PC, 'PC1 falsa cor');
// Map.addLayer(componentesPrincipais.select('pc2'), vis.vis_pc, 'PC2 simple');
param.bandas_PCA.forEach(
    function(momepc){
        print("load componente ", nomepc);
        var temp_pc = componentesPrincipais.select(nomepc);
        Map.addLayer(temp_pc, vis.vis_pc, nomepc + 'simple');
    }
)