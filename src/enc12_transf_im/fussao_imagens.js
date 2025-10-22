/// ========================================================//
// autor : 
// objetivos
// saida: 
//// ======================================================//

/// visualçização 
var vis = {
    // Configurações de visualização
    mosaico_RGB: {
        bands: ['red', 'green', 'blue'],
        min: 0.03,
        max: 0.15,
        gamma: 1.1
    },
    col_RGB: {
        bands: ['B4', 'B3', 'B2'],
        min: 0.03,
        max: 0.15,
        gamma: 1.1
    },
    mosaico_HSI: {
        bands: ['red', 'green', 'blue'],
        min: 0.0,
        max: 0.3,
        gamma: 1.1
    },
    mosaico_brov: {
        bands: ['red', 'green', 'blue'],
        min: 0.0,
        max: 0.025,
        gamma: 1.1
    }
}
/// definir os parametros de entrada 
var param = {
    asset_landsat: "LANDSAT/LC08/C02/T1_RT_TOA",
    // asset_id: "LANDSAT/LC08/C02/T1_RT_TOA/LC08_217071_20230713",
    asset_id: '',
    asset_municipiosBA: "users/CartasSol/shapes/municipios_ba",
    bands_names: ["blue", "green", "red", "nir", "swir1", "siwr2"],
    lst_bnds : ["B2","B3","B4","B5","B6","B7", 'B8'],
};  

// ===========================
// entreda: 
// @multispectral: imagem (6 bandas)
// @panchromatic: imagem (image)
// @saida: image
//==========================================
// Função para aplicar pansharpening usando o método Brovey
function broveyPansharpen(multispectral) {
    // nova_banda_i =  (band_i / ΣBands) × Panchromatic
    // Normalizar as bandas multiespectrais  (band_i / ΣBands) 
    // panchromatic
    var lst_bnd = ["B2","B3","B4","B5","B6","B7"];
    var panchromatic = multispectral.select(['B8'])
    multispectral = multispectral.select(lst_bnd).rename(param.bands_names);
    var raster_sum = multispectral.reduce(ee.Reducer.sum());   // ΣBands
    var normalized = ee.Image(multispectral).divide(raster_sum);  // band_i / ΣBands
    
    // Multiplicar cada banda normalizada pela banda pancromática
    var sharpened = normalized.multiply(panchromatic);

    return sharpened;
}

// Função para aplicar pansharpening usando o método IHS
function ihsPansharpen(multispectral) {
    // Selecionar bandas RGB para transformação IHS
     var lst_bnd = ["B2","B3","B4"];
     var bandsRGB = ['blue', 'green', 'red'];
    var panchromatic = multispectral.select(['B8'])
    multispectral = multispectral.select(lst_bnd).rename(bandsRGB);
    var rgb = multispectral.select(['red', 'green', 'blue']);
    
    // Converter RGB para IHS
    var ihs = rgb.rgbToHsv();
    
    // Substituir o componente Intensity (banda V) pela banda pancromática
    var intensityReplaced = ihs.select(['hue', 'saturation'])
                                .addBands(panchromatic.rename('value'));
    
    // Converter de volta para RGB
    var sharpened = intensityReplaced.hsvToRgb();
    // print("mostrar metadados sharpened IHS ", sharpened);
    
    return sharpened;
}

function export_image(raster_img, name_exp, mygeom){
    var pmto = {
        image: raster_img,
        description: name_exp,
        scale: 15,
        folder: "pansharpening",
        region: mygeom.buffer(5000).bounds(),
        maxPixels: 1e13 
    }
    // Exportar resultado (opcional)
    Export.image.toDrive(pmto);
    print("exportando imagem " + name_exp + " para o folder ")
}

/////////////////==========================================//
// Parte principal 
///======================================================///

//==== carregar Dados
var imagem = null; 
var multi_img = false;
var export_img = false;
var shp_munBahia = ee.FeatureCollection(param.asset_municipiosBA)
                                .map(function(feat){return feat.set('id_cod', 1) });

var mask_ba = shp_munBahia.reduceToImage(['id_cod'], ee.Reducer.first());
if (param.asset_id !== ""){
    imagem = ee.Image(param.asset_id);
}else{
    var data_inicio = ee.Date.fromYMD(2023, 1, 1);
    
    // Carregar imagem Landsat 8
    var img_col_l8 = ee.ImageCollection(param.asset_landsat)
                        .filterBounds(shp_munBahia)
                        .filterDate(data_inicio, data_inicio.advance(1, "year"))
                        .filter(ee.Filter.lt("CLOUD_COVER_LAND", 20))
                        .sort("CLOUD_COVER_LAND", false);
    print("show metadados Img Collection ", img_col_l8.limit(3));
    // print("size imgCollection", img_col_l8.size());
    multi_img = true;
    // imagem = ee.Image(img_col_l8.first());
}

if (multi_img){
    // Extrair bandas multiespectrais (30m) e pancromática (15m)
    var col_multispectral = img_col_l8.select(param.lst_bnds);

    // Aplicar pansharpening ou fussão de imagens 
    var pansharpenedBrovey = col_multispectral.map(broveyPansharpen);
    var pansharpenedIHS = col_multispectral.map(ihsPansharpen);
    print("mostrar metadados de pansharpenedBrovey",  pansharpenedBrovey );
    print("mosatra metadados de  pansharpenedIHS", pansharpenedIHS);
    var multispectral = img_col_l8;
    img_col_l8 = img_col_l8.mosaic().updateMask(mask_ba);
    pansharpenedBrovey = pansharpenedBrovey.mosaic().updateMask(mask_ba);
    pansharpenedIHS = pansharpenedIHS.mosaic().updateMask(mask_ba);
}else{
    print("Mostrar metadado de imagem ", imagem);

    // Extrair bandas multiespectrais (30m) e pancromática (15m)
    var multispectral = imagem.select(param.lst_bnds);

    // Aplicar pansharpening ou fussão de imagens 
    var pansharpenedBrovey = broveyPansharpen(multispectra);
    var pansharpenedIHS = ihsPansharpen(multispectral);
}



// Adicionar camadas ao mapa
if(multi_img){
    Map.addLayer(img_col_l8, vis.col_RGB, 'Original Multiespectral (30m)');
}else{
    Map.addLayer(multispectral, vis.mosaico_RGB, 'Original Multiespectral (30m)');
}
Map.addLayer(pansharpenedBrovey, vis.mosaico_brov, 'Pansharpened Brovey (15m)');
Map.addLayer(pansharpenedIHS, vis.mosaico_HSI, 'Pansharpened IHS (15m)');
// Map.centerObject(imagem.geometry(), 10);

// Comparar resoluções
if (multi_img){
    print('Resolução Original Multiespectral:', multispectral.first().projection().nominalScale());
    print('Resolução Pansharpened:', pansharpenedBrovey.first().projection().nominalScale());
}
if (export_img){
    // WRS_PATH: 216     WRS_ROW:  66
    var name_export = 'L8_path_216_row_66_pansharpened_Brovey_2023';
    export_image(pansharpenedBrovey, name_export, imagem.geometry());
}

