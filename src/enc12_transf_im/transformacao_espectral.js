
// definir parametros
var vis = {
    mosaico: {
        min: 0.03, max: 0.16,
        bands: ["red","green","blue"]
    },
    mos_hsv: {
        min: 0.0, max: 1,
        bands: ["hue","saturation","value"]
    },
    // Defina os parâmetros de visualização para cada componente.
    tasseled_cap: {
        min: -0.1, max: 0.5, 
    },
    ind_veg : {
        min: 10400, max: 18000,
        palette: ["#f09e51", "#f0cb58", "#a8e4a0", "#0e630e"],
    },
    ind_water: {
        min: 2300, max: 10000,
        palette: ["#a8e4a0", "#ffffac", "#6fa8dc", "#3179fd"],
    },
    vis_pc:{
        min: -0.49, max: 2.00,
    },
}
var param = {
    asset_landsat: "LANDSAT/LC08/C02/T1_RT_TOA",
    asset_id: "LANDSAT/LC08/C02/T1_RT_TOA/LC08_217071_20230713",
    asset_municipiosBA: "users/CartasSol/shapes/municipios_ba",
    bands_names: ["blue", "green", "red", "nir", "swir1", "siwr2"],
    lst_bnds : ["B2","B3","B4","B5","B6","B7"],
};  


// ====== function  =======================

function mostrar_histograma(img, lsta_bnd, num_bnd){
    img =  ee.Image(img);
    var lst_cores = [];
    if (num_bnd > 3){
        lst_cores = ['#006b99', '#00af07', '#ef5134', '#6d6b99', '#f0af07', '#c459ae'];
    }else{
        lst_cores = ['#006b99', '#00af07', '#ef5134'];
    }
    // Define the chart and print it to the console.
    var chart = ui.Chart.image.histogram({
        image: img, 
        region: img.geometry(), 
        scale: 360
    })
    .setSeriesNames(lsta_bnd)
    .setOptions({
    title: 'image LANDSAT SR Reflectance Histogram',
    hAxis: {
        title: 'Reflectance (x1e4)',
        titleTextStyle: {italic: false, bold: true},
    },
    vAxis:
        {title: 'Count', titleTextStyle: {italic: false, bold: true}},
    colors: lst_cores
    });
    print("histograma da imagem ", chart);

}

function calcular_tasseled_cap(img, queLandsat){
    // Defina a matriz de coeficientes TCT com base na sua tabela.
    // Vamos usar um array 2D, onde cada linha representa uma componente (Brightness, Greenness, Wetness).
    // Os valores são os coeficientes para cada banda.
    var coeficientesTC = null;

    if (queLandsat !== 8){
        coeficientesTC =  [
            [ 0.3029,  0.2786,  0.4733,  0.5599,  0.5080,  0.1872],  // TC Brightness
            [-0.2941,  -0.2430, -0.5424,  0.7276,  0.0713, -0.1608], // TC Greeness
            [ 0.1511,   0.1973,  0.3283,  0.3407, -0.7117, -0.4559] // TC Wetness
        ];
    }else{
        coeficientesTC= [
            [ 0.3037,  0.2793,  0.4743,  0.5585,  0.5082,  0.1863],  // TC Brightness
            [-0.2848, -0.2435, -0.5436,  0.7243,  0.0840, -0.1800],  // TC Greeness
            [ 0.1509,  0.1973,  0.3279,  0.3406, -0.7112, -0.4572]   // TC Wetness
        ];
    }
    coeficientesTC = ee.Array(coeficientesTC);

    var imageArray = img.toArray().toArray(1);
    print("Show image ", img);
    // Aplique a transformação TCT.
    // A função .reduce() realiza a multiplicação da matriz de bandas pela matriz de coeficientes.
    var bnd_names_tc = ['Brightness', 'Greenness', 'Wetness'];
    // Aplica a transformação TCT.
    var tasseledCap = ee.Image(coeficientesTC).matrixMultiply(imageArray).arrayProject([0]);
    tasseledCap = tasseledCap.arrayFlatten([bnd_names_tc]);
    print("Show Tasseled Cap ", tasseledCap);
    
    return tasseledCap; 
}

function calcular_spectral_index(img){
    // 
    var ind_ndvi = img.normalizedDifference(["nir", "red"]);
    ind_ndvi =  ind_ndvi.add(1).multiply(10000).rename("ndvi");
    // 
    var ind_ndwi = img.normalizedDifference(["green", "nir"]);
    ind_ndwi =  ind_ndwi.add(1).multiply(10000).rename("ndwi");
    // SAVI
    // vSAVI = ((NIR - Red) / (NIR + Red + L)) * (1 + L)
    var ind_savi = img.expression(
        "float( ((b('nir') - b('red')) / (b('nir') + b('red') + 0.5) ) * (1 + 0.5))"
    );
    // var ind_savi = img.expression(
    //     "float( ((nir - red) / (nir + red + L) ) * (1 + L))",
    //     {
    //         nir: "nir",
    //         red: "red",
    //         L: 0.5
    //     }
    // )
    ind_savi = ind_savi.add(1).multiply(10000).rename("savi");

    return img.addBands(ind_ndvi).addBands(ind_ndwi)
                .addBands(ind_savi);
}


// =================== ÍNDICES (NDBSI + NDVI + LST + Wetness) ===================
function building_PCA_space(image, my_region){

    var bandNames = [
                "blue", "green", "red", "nir", "swir1", "siwr2",
                "ndvi","ndwi","savi","Brightness","Greenness", 
                "Wetness"
        ];  


    // Esta função é a que vai retornar a PC em uma nova imagem.
    var getPrincipalComponents = function(centered) {
        
        var namePCA = [
                'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6',
                'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12',
            ];
        var name_sd = [
            'sd1', 'sd2', 'sd3', 'sd4', 'sd5', 'sd6',
            'sd7', 'sd8', 'sd9', 'sd10', 'sd11', 'sd12'        
        ];

        // Colapsa as bandas da imagem em array 1D por pixel.
        var arrays = centered.toArray();
        
        // Calcula a covariancia das bandas na area de estudo.
        var covar = arrays.reduceRegion({
            reducer: ee.Reducer.centeredCovariance(),
            geometry: ee.Geometry(my_region),
            scale: 30,
            maxPixels: 1e13,
            tileScale: 8
        });

        // 'cast' o array de covariancia como array.
        // isto vai representar a covariancia banda a banda dentro da região.
        var covarArray = ee.Array(covar.get('array'));

        // calcula os eigen vectors e eigen values, depois separa os dois.
        var eigens = covarArray.eigen();

        // Eigenvalues.
        var eigenValues = eigens.slice(1, 0, 1);
        
        // Calcula o percentual de variância pra cada componente
        // Assim podemos decidir quantos compnentes vamos manter
        var eigenValuesList = eigenValues.toList().flatten()
        var total = eigenValuesList.reduce(ee.Reducer.sum())

        var percentageVariance = eigenValuesList.map(
                function(item) {
                    var component = eigenValuesList.indexOf(item).add(1).format('%02d')
                    var variance = ee.Number(item).divide(total).multiply(100).format('%.2f')
                    return ee.List([component, variance])
            }
        )
    
        // Cria um dicionário que usamos pra setar
        // as propriedades na imagem final
        var varianceDict = ee.Dictionary(percentageVariance.flatten())
        
        // Cria a matriz dos eigenvectors como linhas.
        var eigenVectors = eigens.slice(1, 1);
        
        // Converte o array original para um array 2D para usar no calculo de matrizes.
        var arrayImage = arrays.toArray(1);

        // Multiplica a imagem array pela matriz de eigenvectors.
        var principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage);

        // Transforma a raiz quadrada dos eigenvalues em uma imagem.
        // Usa o abs() pra converter valores negativos pra positivos
        // Acha a raiz quadrada
        var sdImage = ee.Image(eigenValues.abs().sqrt())
                            .arrayProject([0])
                                .arrayFlatten([name_sd]);

        // Converte as PCs em uma imagem normalizada pelo Desvio Padrão.
        return principalComponents
                    // Joga fora qualquer dimensão indesejada, [[]] -> [].
                    .arrayProject([0])
                    // Transforma o array de uma banda em uma multi-band image, [] -> image.
                    .arrayFlatten([namePCA])
                    // Normaliza pelos desvios padrão.
                    .divide(sdImage)
                    .set(varianceDict);
    };

    
    
    // Calcula a média dos dados (acelerar a redução da matriz de covariancia)
    var meanDict = image.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: ee.Geometry(my_region),
        scale: 30,
        maxPixels: 1e13,
        tileScale: 8
    });
    
    // Subtrai cada banda do seu valor médio
    var means = ee.Image.constant(meanDict.values(bandNames));
    var centered = image.subtract(means);
    
    // make image array PCA 
    var pcImage = getPrincipalComponents(centered);
    pcImage = pcImage.select(['PC1', 'PC2', 'PC3']);
    return image.addBands(pcImage);
}


// Function to add an NDVI band, the dependent variable.
function agregateBandsgetFractions (image) {
    // Define endmembers
    var endmembers =  [
        [0.05, 0.09, 0.04, 0.61, 0.30, 0.10], //*gv*/
        [0.14, 0.17, 0.22, 0.30, 0.55, 0.30], //*npv*/
        [0.20, 0.30, 0.34, 0.58, 0.60, 0.58], //*soil*/
        [0.0 , 0.0,  0.0 , 0.0 , 0.0 , 0.0 ], //*Shade*/
        [0.90, 0.96, 0.80, 0.78, 0.72, 0.65] //*cloud*/        
    ]

    // Uminxing data
    var bandas = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']
    var bandsFraction = ['gv','npv','soil','shade','cloud']

    var fractions = ee.Image(image).select(bandas)
                    .unmix({
                        endmembers: endmembers
                    }).float()
    
    fractions = fractions.select([0,1,2,3,4], bandsFraction)
    // GVshade = GV /(1 - SHADE)
    // NDFIa = (GVshade - SOIL) / (GVshade + )
    // ((GV / (1 - SHADE)) - SOIL) / ((GV / (1 - SHADE)) + (NPV) + SOIL)
    var ndfia = fractions.expression(
        "float((b('gv') / (1 - b('shade'))) - b('soil')) / float((b('gv') / (1 - b('shade'))) + b('npv') + b('soil'))")
    
    ndfia = ndfia.toFloat().rename('ndfia')    
    return image.addBands(fractions).addBands(ndfia)
};

//==== carregar Dados
var imagem = null; 
var show_plot = false;
var calcularPCA = true;
var shp_munBahia = ee.FeatureCollection(param.asset_municipiosBA);
if (param.asset_id !== ""){
    imagem = ee.Image(param.asset_id);
}else{
    var data_inicio = ee.Date.fromYMD(2023, 1, 1);
    var img_col_l8 = ee.ImageCollection(param.asset_landsat)
                            .filterDate(data_inicio, data_inicio.advance(1, "year"))
                            .filterBounds(pto1)
                            .filter(ee.Filter.lt("CLOUD_COVER_LAND", 30))
                            .sort("CLOUD_COVER_LAND" );
    print("show metadados Img Collection ", img_col_l8.limit(3));
    print("size imgCollection", img_col_l8.size());

    imagem = ee.Image(img_col_l8.first())
}
        
// seleccionar as 6 primeiras bandas e renomear elas
imagem = imagem.select(param.lst_bnds).rename(param.bands_names);
print("show metadata Image ", imagem);

/// transformacao HSV
var bnd_rgb = ["blue", "green", "red"];
var imagem_hsv = imagem.select(bnd_rgb).rgbToHsv();
print("mostrar imagem ", imagem_hsv);

if (show_plot){
    print("mostrar histrograma das bandas RGB,Nir,s1,s2");
    mostrar_histograma(imagem, param.bands_names, 6);
    print("mostra histograma da imagem HSV");
    mostrar_histograma(imagem_hsv, ["hue","saturation","value"], 3);
}

var mTasseledCap = calcular_tasseled_cap(imagem, 8);

//  calcular indices espectrais 
imagem = calcular_spectral_index(imagem);
// juntando todas as bandas 
imagem = imagem.addBands(mTasseledCap);
print("ver lista de bandas ", imagem);

var raster_pca = null;
if (calcularPCA){
    raster_pca = building_PCA_space(imagem, imagem.geometry())
    print("show metadado do raster PC", raster_pca);
    raster_pca = raster_pca.select(['PC1', 'PC2', 'PC3'])
}


Map.addLayer(shp_munBahia, {color: "green"}, "miunicipiosBA", false);
Map.addLayer(imagem, vis.mosaico, "imagem", true);
Map.addLayer(imagem_hsv, vis.mos_hsv, "imagem hsv", false);
Map.addLayer(mTasseledCap.select('Brightness'), vis.tasseled_cap, 'TC - Brightness', false);
Map.addLayer(mTasseledCap.select('Greenness'), vis.tasseled_cap, 'TC- Greenness', false);
Map.addLayer(mTasseledCap.select('Wetness'), vis.tasseled_cap, 'TC - Wetness', false);
Map.addLayer(imagem.select("ndvi"), vis.ind_veg, "raster NDVI", true);
Map.addLayer(imagem.select("ndwi"), vis.ind_water, "raster NDWI", true);
Map.addLayer(imagem.select("savi"), vis.ind_veg, "raster SAVI", true);


if (calcularPCA){
    var lst_bdnPC = ['PC1', 'PC2', 'PC3'];
    lst_bdnPC.forEach(function(bnd_name){
        print("show band ==>  " + bnd_name);
        Map.addLayer(raster_pca.select(bnd_name), vis.vis_pc, bnd_name, true);
    })

}