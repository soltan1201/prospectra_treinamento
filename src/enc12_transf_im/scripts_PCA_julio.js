
var municipios = ee.FeatureCollection('projects/ee-my-jpedrassoli/assets/RECORTES-GEOGRAFICOS/BR_Municipios_2021_com_regioes')
var feiraSantana = municipios.filter(ee.Filter.eq('NM_MUN', 'Feira de Santana'))


function PCA(maskedImage, area){
    var image = maskedImage.unmask()
    var bandNames = image.bandNames()
    
    // Calcula a média dos dados (acelerar a redução da matriz de covariancia)
    var meanDict = image.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: area,
        scale: 30,
        maxPixels: 1e13,
        tileScale: 16
    });
    
    // Subtrai cada banda do seu valor médio
    var means = ee.Image.constant(meanDict.values(bandNames));
    var centered = image.subtract(means);
    
    // Retorna os nomes das bandas novas
    var getNewBandNames = function(prefix) {
        var seq = ee.List.sequence(1, bandNames.length());
        return seq.map(
                function(b) {
                    return ee.String(prefix).cat(ee.Number(b).int());
                });
    };
    
    // Esta função é a que vai retornar a PC em uma nova imagem.
    var getPrincipalComponents = function(centered) {
        
        // Colapsa as bandas da imagem em array 1D por pixel.
        var arrays = centered.toArray();
        
        // Calcula a covariancia das bandas na area de estudo.
        var covar = arrays.reduceRegion({
                                reducer: ee.Reducer.centeredCovariance(),
                                geometry: area,
                                scale: 30,
                                maxPixels: 1e13,
                                tileScale: 16
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
                                .arrayProject([0]).arrayFlatten([getNewBandNames('sd')]);

        // Converte as PCs em uma imagem normalizada pelo Desvio Padrão.
        return principalComponents
                    // Joga fora qualquer dimensão indesejada, [[]] -> [].
                    .arrayProject([0])
                    // Transforma o array de uma banda em uma multi-band image, [] -> image.
                    .arrayFlatten([getNewBandNames('pc')])
                    // Normaliza pelos desvios padrão.
                    .divide(sdImage)
                    .set(varianceDict);
    };
    
    var pcImage = getPrincipalComponents(centered);
    
    return pcImage.mask(maskedImage.mask());
}
