# AnalyseSpecimenColours
# This repository contains:

## (1) R Package to execute colourimetric analyses as described in Teo et al. (2025)

### (2) Supplementary materials of raw colour data for Teo et al. (2025)

#### To install (1) in R

	install_github("dbssyck/AnalyseSpecimenColours")

### If you use this package please cite
- Teo, S.H., Sin, Y.C.K., Nieves, M.M.I., Rheindt, F.E., 2025. Birds of a feather: Comprehensive plumage colour analysis for a revised subspecies classification of the Chestnut-winged Babbler (Cyanoderma erythropterum) species complex. Avian Research, 100307. https://doi.org/10.1016/j.avrs.2025.100307

This script uses the 75% subspecies rule described in Amadon (1949) to quantify whether two populations are sufficiently differentiated to warrant being classified as separate subspecies. You can modify the thresholds to apply the 95% (Patten & Unitt, 2002) subspecies rule as well. The following two works should be cited accordingly.
- Amadon, D., 1949. The seventy-five per cent rule for subspecies. Condor, 51(6), 250–258. https://doi.org/10.2307/1364805 
- Patten, M.A., Unitt, P., 2002. Diagnosability versus mean differences of Sage Sparrow subspecies. The Auk, 119(1), 26-35. https://doi.org/10.1093/auk/119.1.26

The script also uses the CIEDE 2000 formula (Alman et al., 2001) to quantify colour differences between specimens.
- Alman, D.H., Komatsubara, H., Li, W., Luo, M.R., Melgosa, M., Nobbs, J.H., Rigg, B., Robertson, A. R., Witt, K., Cui, G., 2001. Improvement to Industrial Colour-Difference Evaluation. (CIE Pub. No. 142.) International Commission on Illumination, Vienna, Austria. ISBN: 978 3 901906 08 4

Finally, this script relies on the following external packages that should be acknowledged accordingly
- colorscience
- MASS
- ape
- ggplot2
- openxlsx

## About

The 75% subspecies rule proposed by Amadon (1949) states that in order for two populations to be considered distinct subspecies, for any given taxonomically important characteristic, the chracter measurements for 75% of the individuals from a population must not overlap with the character measurements of X% of the individuals from the other population. The specific percentage X% can range from 75% to 99%. A more rigorous threshold has also been proposed by Patten and Unitt (2002), whereby 95% of a population must not overlap with 99% of the other. Museum specimens serve as crucial data source that allow us to conduct such taxonomic assessments.

While the overlap (or lack of thereof) in characteristics such as beak and tarsus length can be easily computed, overlaps in colour differences are much more difficult to interpret. Colour spectrometry has traditionally been the standard for measuring the colour differences in specimen. In Teo et al. (2025) we show that such differences can be calculated from photographs as well if the specimens are all photographed in a consistent manner (specimen positioning, lighting, camera settings etc).

From consistently taken photographs, the average RGB values of taxonomically meaningful body parts can be measured (using external programs such as ImageJ or MATO).

#### 1) This package uses the RGB values to test the 75% subspecies rule between specimens belonging to different populations by first conducting a Linear Discriminant Analysis (LDA) on the RGB values, then running the Amadon (1949) subspecies rule on the LD1 values.

#### 2) This package also compute the colour differences between each pair of specimens using the CIEDE2000 colour distance metric, then generate a phylogram that can aid in determining how the clustered the specimens are in terms of their characters


## Usage

To use this package simply run:

    run.functions(colour.data, accession.number, grouping.column,
                  R.col.name, G.col.name, B.col.name,
                  run = "key",
                  threshold.1 = 0.75, threshold.2 = 0.99,
                  export.colours.plot = FALSE, export.figures = TRUE, export.excel = TRUE,
                  output.folder = NULL, excel.file.name = NULL)

### Arguments
##### colour.data
Data frame containing sample information and RGB values of each sample. This data frame can have however many columns or rows. Importantly, each row must correspond to a sample, and there must be 1) at least one column providing each sample a unique accession number (unique dummy values can be used if the specimens did not have accession numbers), and 2) at least one column with the population that the samples are assigned to (e.g. subspecies, geographic locality). See example data frame below.

##### accession.number
Name of column with accession number

##### grouping.column
Name of column with population assignment. You can choose to group your samples in different ways depending on your scope. In Teo et al. (2025), we ran separate analyses for "Subspecies" and "General Locality". The latter was conducted to explore if there were any differences between specimens currently considered the same subspecies despite coming from different locations.

##### R.col.name, G.col.name, B.col.name
Names of columns containing the R values, G values, and B values respectively

##### run
"key" or "all". "key" runs the key analyses found to be important in Teo et al. (2025), namely running the 75% subspecies rule test between each pair of populations on the LD1 values after conducting LDA on the RGB values, and also generating a phylogram based on CIEDE2000 colour differences between each pair of samples. "all" runs all other exploratory analyses described in Teo et al. (2025), which include running the same set of analyses after converting the RGB values to different colour models (CMY, HSV, and HSL), and also running the 75% subspecies rule on the individual colour values (such as R, G, and B), principal components, and linear discriminants of each colour model.

##### threshold.1
Default value 0.75. Threshold value of the percentage of a population that must be separable from threshold.2 of the other population. To test the Patten and Unitt (2002) rule, change this value to 0.95 and keep threshold.2 as 0.99.

##### threshold.2
Default value 0.99. Threshold value of the percentage of a population that must be separable from threshold.1 of the other population.

##### export.colours.plot
TRUE/FALSE. Defaults to FALSE. TRUE exports boxplots of individual colour values (e.g. R, G and B) for all colour models, for all pairwise group comparisons. However these plots will only be generated if run = "all" as well.

##### export.figures
TRUE/FALSE. Defaults to FALSE. TRUE exports the LDA and PCA figures into an output file.

##### export.excel
TRUE/FALSE. TRUE exports results of the subspecies rule test as an .xlsx file.

##### output.folder
Path to folder where the image and excel files will be exported. If not designated, all files will be written in the present working directory. If a folder with the name designated here doesn't exist in the present working directory, it will be created.

##### excel.file.name
Specify file name of excel file (if export.excel = TRUE). If undefined, the file will be named diagnosability.xlsx

### Explanation of output excel sheet
##### all.results.summary
Summarises whether each pair of populations are diagnosable from one another using the subspecies rule test. Each column corresponds to the characters that were tested.

All the other sheets contain specific results for each tested character, detailing 1) which population within the pair has a lower mean value, 2) whether the lower population is diagnosable from the higher population (i.e. if the top 25% of the individuals in the lower population only overlaps with the bottom 1% of individuals in the higher population or less, the result will be "Yes", and anything beyond will be considered undiagnosable = "No"), 3) whether the higher population is diagnosable from the lower population (i.e. if the bottom 25% of the individuals in the higher population only overlaps with the top 1% of individuals in the lower population or less, the result will be "Yes", and anything beyond will be considered undiagnosable = "No"), and 4) whether they are diagnosable from each other (i.e. both lower.from.higher and higher.from.lower are "Yes").

## If you want to explore your dataset in a step wise fashion...

run.functions() runs the full set of analyses described in Teo et al. (2025). The functions that were written for this package and are called in run.functions() can be ran individually as well.

### Pairwise Delta
Calculates the CIEDE2000 values between every pair of samples in the dataset.

    pairwise.delta(colour.data, accession.number, R.col, G.col, B.col, as.distance.matrix)

##### colour.data, accession.number, R.col, G.col, B.col
As above

##### as.distance.matrix
TRUE/FALSE. Defaults to FALSE. TRUE outputs the pairwise distance in R's distance matrix format. FALSE outputs pairwise distances as a mirrored matrix.

### Testing subspecies rule

Runs the subspecies rule test between two populations. You can use this function if testing the subspecies rule between more straightforward measurements such as wing length, bill depth etc.

    check.linear.overlap.Amadon(data.1, data.2, name.data.1, name.data.2, threshold.1, threshold.2)

##### data.1
Values from first population
##### data.2
Values from second population
##### name.data.1
Name of first population
##### name.data.2
Name of second population
##### threshold.1, threshold.2
As above

### Running subspecies rule test between population pairs

Runs check.linear.overlap.Amadon() on all pairwise combinations of populations.

    pairwise.group.Amadon(colour.data, accession.number, population, parameter.to.check.overlap, threshold.1, threshold.2)

##### colour.data
As above. If running this function for characters such as wing length, those data should be included in this data frame rather than having RGB values. 
##### accession.number, population
As above
##### parameter.to.check.overlap
Name of column with character values to test the subspecies rule
##### threshold.1, threshold.2
As above

### Run neighbour-joining trees from pairwise CIEDE2000 values

Creates neighbour-joining tree based on RGB input and produce a phylogram

    neighbour.joining(colour.data, accession.number, population, R.col, G.col, B.col, export.figures, output.folder, output.file.prefix)

##### colour.data, accession.number, population, R.col, G.col, B.col
As above
##### export.figures
TRUE/FALSE. Defaults to FALSE. TRUE exports the phylogram as an image file.
##### output.folder
As above
##### output.file.prefix
Specify file name of phylogram image file (if export.figures = TRUE). If undefined, the file will be named nj_phylogram.jpg

### Convert to other colour models

Convert RGB colour model to greyscale, cmy, hsl, and hsv using the colorscience package

    add.colour.models.from.RGB(colour.data, R.col, G.col, B.col)

##### colour.data, R.col, G.col, B.col
As above

### Run subspecies rule test after conducting PCA and LDA on colour values

Test subspecies rule after running PCA and LDA on the colour values

    pca.lda.overlap(colour.data, accession.number, population, R.col, G.col, B.col, threshold.1, threshold.2, colour.data.from.all.colour.models, export.figures, output.folder) 

##### colour.data, accession.number, population, R.col, G.col, B.col, threshold.1, threshold.2
As above
##### colour.data.from.all.colour.models
If undefined, the subspecies rule test will be conducted on RGB values only. If testing the subspecies rule for other colour models as well, put the outputs from add.colour.models.from.RGB() here
export.figures, output.folder: As above

## Example of data format

| Specimen   ID | General Locality |             Specific Locality            | Year |   Subspecies  |  R  |  G  |  B |
|:-------------:|:----------------:|:----------------------------------------:|:----:|:-------------:|:---:|:---:|:--:|
| ZRC.3.21483   | Sumatra          | Sungei Tasik, Langkat, Sumatra           | 1919 | pyrrhophaeum  |  77 |  65 | 49 |
| ZRC.3.21442   | S. Pen. Malaysia | Kuala Tembeling,   Pahang                | 1907 | erythropterum | 116 | 102 | 78 |
| ZRC.3.21473   | Natuna           | Bungurun, North   Natuna                 | 1928 | neocarum      |  93 |  80 | 60 |
| ZRC.3.21457   | Tioman           | Juara Bay, Pulau   Tioman                | 1915 | erythropterum |  80 |  66 | 50 |
| ZRC.3.21402   | Sabah            | Rayoh, North Borneo                      | 1928 | bicolor       |  87 |  74 | 54 |
| ZRC.3.21371   | Sabah            | Semawang Road,   Sandakan                | 1927 | bicolor       |  76 |  66 | 52 |
| ZRC.3.21381   | Banggi           | Banguay Island, North   Borneo           | 1927 | bicolor       |  94 |  80 | 61 |
| ZRC.3.21406   | Sarawak          | Lingit, Saribas,   Sarawak               | 1919 | bicolor       |  79 |  69 | 55 |
| ZRC.3.21481   | Natuna           | Bungurun, North   Natuna                 | 1928 | neocarum      |  76 |  65 | 51 |
| ZRC.3.21423   | Thailand         | Ban Kok Klap, Bandon,   Malay Peninsular | 1913 | sordidum      |  86 |  75 | 58 |
| ZRC.3.21491   | Sumatra          | Tanjong Kanan,   Sumatra                 | 1939 | pyrrhophaeum  |  81 |  67 | 53 |
| ZRC.3.21407   | Sarawak          | Saribas, Sarwak                          | 1922 | bicolor       |  69 |  55 | 43 |
| ZRC.3.21493   | Sumatra          | Tanjong Kanan,   Sumatra                 | 1939 | pyrrhophaeum  |  82 |  68 | 54 |
| ZRC.3.21385   | Sabah            | Kudat, North Borneo                      | 1927 | bicolor       |  80 |  70 | 56 |
| ZRC.3.21379   | Sabah            | Bettotan, Sandakan,   Borneo             | 1927 | bicolor       |  86 |  75 | 61 |
| ZRC.3.21396   | Sabah            | Rayoh, North Borneo                      | 1928 | bicolor       |  74 |  64 | 52 |
| ZRC.3.21433   | N. Pen. Malaysia | Tanjong Antu,   Dindings                 | 1918 | erythropterum |  95 |  82 | 60 |
| ZRC.3.21422   | Thailand         | Chao Ram, Siam                           | 1922 | sordidum      |  99 |  85 | 63 |
| USNM 174799   | Natuna           | Pulau Natuna Besar                       | 1900 | neocara       | 102 |  85 | 77 |
| ZRC.3.21466   | S. Pen. Malaysia | Kampung Kudu,   Southwest Johor          | 1904 | erythropterum | 123 | 113 | 85 |
