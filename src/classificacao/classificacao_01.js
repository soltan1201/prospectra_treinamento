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
    bands: ["B4", "B3", "B2"],
  },
  classificacao: {
    min: 0,
    max: 5,
    palette: paletaClasses,
  },
};

/*
Etapa 01: Seleção dos municípios e definição dos contornos

Seleciona apenas o estado e depois realiza um novo filtro para cidades dentro desse estado
Contorno é definido com style, para não ficar uma cor sólida
Variável municipiosSelecionados é gerenciada pelo gee para que as três imagens desconexas se tornem
uma feature collection
*/

var lista_mun = ["Barreiras", "Remanso", "Jacobina"];

var municipiosSelecionados = table
  .filter(ee.Filter.eq("ADM1_NAME", "Bahia"))
  .filter(ee.Filter.inList("ADM2_NAME", lista_mun));

// Adiciona o contorno de cada município separadamente
lista_mun.forEach(function (nomeMunicipio) {
  var munIndividual = municipiosSelecionados.filter(
    ee.Filter.eq("ADM2_NAME", nomeMunicipio)
  );
  Map.addLayer(munIndividual.style(vis.mun), {}, "Município: " + nomeMunicipio);
});

print(municipiosSelecionados);

/*
Etapa 02: Composite

Criar uma única imagem
Geometry definido faz com que o processamento se limite somente as áreas escolhidas
Filtros:
Bounds > Descarta imagens do satélite que não passaram pelos meus polígonos
Date > Definição do período de busca
lt cloudy pixel percentage > Descarta imagens com mais de 20% de nuvens
Median > Pega a pilha de imagens que restaram e calcula a mediana de cada pixel, formando uma imagem limpa
clip > corta a imagem exatamente no formato do nosso município
*/

var MAX_PROB_NUVEM = 30;

function mascararNuvens(img) {
  img = ee.Image(img);
  var nuvens = ee.Image(img.get("cloud_mask")).select("probability");
  var naoNuvem = nuvens.lt(MAX_PROB_NUVEM);
  return img.updateMask(naoNuvem);
}

var geometry = municipiosSelecionados.geometry();
Map.centerObject(geometry, 6);

var s2 = ee
  .ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(geometry)
  .filterDate("2023-01-01", "2024-01-01")
  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 30));

var s2Clouds = ee
  .ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
  .filterBounds(geometry)
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
//.clip(geometry);

var mosaico = ee.ImageCollection(filtroS2).median();

print("Mosaico: ", mosaico);

Map.addLayer(mosaico.clip(geometry), vis.mosaico, "Mosaico Sentinel-2 (2023)");

/*
Etapa 03: Coleta de amostras e unificação das classes

CLASSES =>

Água (ID: 0)
Urbano (ID: 1)
Solo Exposto (ID: 2)
Vegetação (ID: 3)
Pastagem (ID: 4)
Agricultura (ID: 5)

Usa função desenhar forma para capturar amostras

Todas as classes são unidas em uma única feature collection
Valores de reflectância são extraídos
Bandas de cada classe são extraídas
*/

var amostrasTreinamento = agua
  .merge(urbano)
  .merge(soloExposto)
  .merge(vegetacao)
  .merge(pastagem)
  .merge(agricultura);

var dadosTreinamento = mosaico.sampleRegions({
  collection: amostrasTreinamento,
  properties: ["classe"],
  scale: 30,
  geometries: false, // CRUCIAL: Não salva o desenho, só os dados
  tileScale: 16, // CRUCIAL: Divide o processamento em 16x para não estourar memória
});

print("Dados de treinamento (size): ", dadosTreinamento.size());

dadosTreinamento = dadosTreinamento
  .randomColumn()
  .filter(ee.Filter.lt("random", 0.2));

//dadosTreinamento = dadosTreinamento.randomColumn().sort('random').limit(5000);

//print('Amostra dos dados:', dadosTreinamento.limit(5));
print("Total de pixels coletados:", dadosTreinamento.size());

// Conta quantos pontos existem para cada valor da propriedade 'classe'
var contagemPorClasse = dadosTreinamento.aggregate_histogram("classe");

print("Quantidade de pixels por classe:", contagemPorClasse);

/*
Etapa 04: Treinamento

Algoritmo escolhido é o RandomForest
Modelo analisa relação entre bandas e classes
Imagem é classificada com paleta de cores
*/

var bandasDeTreinamento = [
  "B2",
  "B3",
  "B4",
  "B5",
  "B6",
  "B7",
  "B8",
  "B11",
  "B12",
];

var classificador = ee.Classifier.smileRandomForest({
  numberOfTrees: 100,
  minLeafPopulation: 1,
});

var classificadorTreinado = classificador.train({
  features: dadosTreinamento,
  classProperty: "classe",
  //inputProperties: mosaico.bandNames()
  inputProperties: bandasDeTreinamento,
});

var imagemClassificada = mosaico.classify(classificadorTreinado);

Map.addLayer(
  imagemClassificada.clip(geometry),
  vis.classificacao,
  "Classificação (Random Forest)"
);

// ANALISAR A IMPORTÂNCIA DAS VARIÁVEIS (CURIOSIDADE TÉCNICA)
// O Random Forest pode nos dizer quais bandas foram mais úteis para diferenciar as classes.
// Isso aparece no Console como um dicionário/gráfico.
var importancia = classificadorTreinado.explain();
print("Importância das Bandas:", importancia);

/*
Etapa 05: Validação

Utilizando índices espectrais para criar um relacionamento com a classificação realizada
*/

function calIndices(img) {
  var ndvi = img.normalizedDifference(["B8", "B4"]).rename("NDVI");
  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  var ndbi = img.normalizedDifference(["B11", "B8"]).rename("NDBI");

  return img.addBands([ndvi, ndwi, ndbi]);
}

var mosaicoIndices = calIndices(mosaico);
print("Mosaico com índices: ", mosaicoIndices);
