//Autor(es): @felipe-py
//Link gee: https://code.earthengine.google.com/fc14c5c468fcc99333a7fbd5f3d9da47

var paletaClasses = [
  "blue", // 0 água
  "red", // 1 urbano
  "yellow", // 2 solo exposto
  "006400", // 3 vegetação (Verde escuro)
  "90EE90", // 4 pastagem (Verde claro)
  "orange", // 5 agricultura
];

var vis = {
  mun: {
    color: "red",
    fillColor: "000000",
    width: 2,
  },
  mosaico: {
    min: 0,
    max: 3000,
    bands: ["red", "green", "blue"],
  },
  classificacao: {
    min: 0,
    max: 5,
    palette: paletaClasses,
  },
};

// Calculos de indices espectrais //
function calculo_ndvi(img) {
  // nome
  var indice_ndvi = img.normalizedDifference(["nir", "red"]).rename("ndvi");
  // nome
  var indice_ndwi = img.normalizedDifference(["green", "nir"]).rename("ndwi");
  // nome
  var indice_ndbi = img.normalizedDifference(["swir1", "nir"]).rename("ndbi");

  return img.addBands([indice_ndvi, indice_ndwi, indice_ndbi]);
}

function export_csv(feat_roi, nome_exp) {
  var pmto_exp = {
    collection: feat_roi,
    description: nome_exp,
    folder: "samples_rois",
  };

  Export.table.toDrive(pmto_exp);
  print("exportando tabela " + nome_exp);
}

var amostrasTreinamento = agua
  .merge(urbano)
  .merge(soloExposto)
  .merge(vegetacao)
  .merge(pastagem)
  .merge(agricultura);

var bnd_S2 = ["B2", "B3", "B4", "B8", "B11", "B12"];
var bnd_names = ["blue", "green", "red", "nir", "swir1", "swir2"];

var lista_mun = ["Barreiras", "Remanso", "Jacobina"];
var asset_municipios = "FAO/GAUL_SIMPLIFIED_500m/2015/level2";
var municipiosSelecionados = ee
  .FeatureCollection(asset_municipios)
  .filter(ee.Filter.eq("ADM1_NAME", "Bahia"))
  .filter(ee.Filter.inList("ADM2_NAME", lista_mun))
  .map(function (feat) {
    return feat.set("idCod", 1);
  });

var MAX_PROB_NUVEM = 30;

function mascararNuvens(img) {
  img = ee.Image(img);
  var nuvens = ee.Image(img.get("cloud_mask")).select("probability");
  var naoNuvem = nuvens.lt(MAX_PROB_NUVEM);
  return img.updateMask(naoNuvem);
}

var s2 = ee
  .ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(municipiosSelecionados.geometry())
  .filterDate("2023-01-01", "2024-01-01")
  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 30));

var s2Clouds = ee
  .ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
  .filterBounds(municipiosSelecionados.geometry())
  .filterDate("2023-01-01", "2024-01-01");

var s2MascaraNuvem = ee.Join.saveFirst("cloud_mask").apply({
  primary: s2,
  secondary: s2Clouds,
  condition: ee.Filter.equals({
    leftField: "system:index",
    rightField: "system:index",
  }),
});

var filtroS2 = s2MascaraNuvem
  .filter(ee.Filter.notNull(["cloud_mask"]))
  .map(mascararNuvens);

var mosaico = ee.ImageCollection(filtroS2).map(function (img) {
  // renomear as bandas da imagem
  return img.select(bnd_S2, bnd_names);
});
print("Mosaico: ", mosaico);
mosaico = mosaico.map(calculo_ndvi);

print(mosaico.first().bandNames());

var bnd_import = ["red", "nir", "swir1", "ndvi"];

// Adiciona o contorno de cada município separadamente
lista_mun.forEach(function (nomeMunicipio) {
  var munIndividual = municipiosSelecionados.filter(
    ee.Filter.eq("ADM2_NAME", nomeMunicipio),
  );
  var mask_mcpio = munIndividual.reduceToImage(["idCod"], ee.Reducer.first());
  var amostrasTreinamentotmp = amostrasTreinamento.filterBounds(
    munIndividual.geometry(),
  );

  var mosaico_mcpio = mosaico
    .filterBounds(munIndividual.geometry())
    .median()
    .updateMask(mask_mcpio);

  var dadosTreinamento = mosaico_mcpio.sampleRegions({
    collection: amostrasTreinamentotmp,
    properties: ["classe"],
    scale: 30,
    geometries: false,
    tileScale: 4,
  });

  print("Dados de treinamento (size): ", dadosTreinamento.size());

  //export_csv(dadosTreinamento, nomeMunicipio + '_rois');

  var amostrasComRandom = dadosTreinamento.randomColumn("random");
  var split = 0.7;
  var trainingPartition = amostrasComRandom.filter(
    ee.Filter.lt("random", split),
  );
  var validationPartition = amostrasComRandom.filter(
    ee.Filter.gte("random", split),
  );

  var arvoresTeste = [50, 100, 300];

  print("--- " + nomeMunicipio + " ---");

  arvoresTeste.forEach(function (nArvores) {
    var modeloTeste = ee.Classifier.smileRandomForest({
      numberOfTrees: nArvores,
      minLeafPopulation: 5,
      seed: 42,
    }).train(trainingPartition, "classe", bnd_import);

    var validacao = validationPartition.classify(modeloTeste);
    var matrix = validacao.errorMatrix("classe", "classification");
    print("Trees: " + nArvores + " Acc:", matrix.accuracy());
  });

  var pmtro_cc = {
    numberOfTrees: 300,
    minLeafPopulation: 5,
    seed: 42,
  };

  // Processo de classificação por municipio
  var classificador = ee.Classifier.smileRandomForest(pmtro_cc).train(
    dadosTreinamento,
    "classe",
    bnd_import,
  );

  var imagemClassificada = mosaico_mcpio.classify(classificador);

  var imagemFinal = imagemClassificada.focalMode({
    radius: 1,
    kernelType: "square",
    units: "pixels",
  });

  Map.addLayer(
    munIndividual.style(vis.mun),
    {},
    "Município: " + nomeMunicipio,
    false,
  );
  Map.addLayer(
    mosaico_mcpio,
    vis.mosaico,
    "Mosaico 2023 (" + nomeMunicipio + ")",
  );
  Map.addLayer(
    imagemFinal,
    vis.classificacao,
    "Classificação (" + nomeMunicipio + ")",
  );

  // ANALISAR A IMPORTÂNCIA DAS VARIÁVEIS (CURIOSIDADE TÉCNICA)
  // O Random Forest pode nos dizer quais bandas foram mais úteis para diferenciar as classes.
  // Isso aparece no Console como um dicionário/gráfico.
  var importancia = classificador.explain();
  print("Importância das Bandas:", importancia);

  var idAsset = "users/luislipecunha/Classificacao_" + nomeMunicipio + "_2023";

  Export.image.toAsset({
    image: imagemFinal,
    description: "Salvar_" + nomeMunicipio,
    assetId: idAsset,
    region: munIndividual.geometry(),
    scale: 10,
    maxPixels: 1e13,
    pyramidingPolicy: { ".default": "mode" },
  });
});
