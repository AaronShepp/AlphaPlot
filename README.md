# AlphaPlot
Layout Optimizer for Designing Species-Specific Pairwise Effect Experiments

## Background 

Species-specific pairwise effect refers to the effect that a neighbor tree has on a focal tree due to their specific species-pair permutation. Species "A" can have a different effect on "B" than "B" has on "A". These relationships can be used to advise species selection for diverse tree planting and to optimize the layout of selected mixtures. There are many important considerations when establishing a common-garden experiment to analyze tree diversity (see Bruelheide et al., 2014 and Vanclay, 2006) This layout optimizer program implements two of my suggestions that prioritize the quantification of species-specific pairwise effects. Firstly, trees should be planted in a hexagonal grid to standardize the distance between neighboring trees. Second, there should be a standardized gradient of neighbor counts for each species pair permutation. In other words, there are at least a certain number of species "A" that are not beside species "B", at least a certain number of "A" that are beside exactly 1 "B", at least a certain number of "A" that are beside exactly 2 "B", etc. These "representation criteria" are applied to every species-pair permutation. The neighbor representation gradient ensures that resulting pairwise effect coefficients describe an effect slope as opposed to a binary effect. Furthermore, a good gradient will enable the creation of a pairwise effect function including the identification of effect saturation that could be highly relevant to forestry. Such a function has never been described, to the best of my knowledge. Finally, by standardizing the neighbor representation gradient you raise the likelihood of successfully quantifying pairwise effects for some species-pairs that would otherwise be under represented in a random layout. This lowers the cost of establishing these experiments by requiring fewer trees to be planted per degree of freedom.

## Requirements
### Python Version
Python 3.8 or later recommended

### Required Packages
pandas — tabular output and Excel export

matplotlib — visualization of planting layouts

The following modules are part of the Python standard library and do not require separate installation: random, os, math, time

### Installing dependencies
If you do not already have the required packages, use the following code in the command line (example for pandas):

```pip install pandas matplotlib```

## Basic Usage
This program is designed to be run as a standalone script. Users can customize experimental parameters. These are defined at the top of the file in the INPUTS section.

### Define plot size
This program creates a hexagonal planting layout where each plant is hexagonally arranged to have 6 neighbors, and the entire plot is an even-sided hexagon. You can define the size of an inner "core" section of trees to be sampled, and an outer "edge" section of trees that serve to minimize edge effects. Core trees will be used as focal trees and neighbor trees for optimizing neighbor representation counts. Edge trees will never be used as focal trees but will be used as neighbor trees if they are adjacent to the core. 

In this README, I have written every line of code where the user can customize parameters. For example, ```INNER_DIAMETER = 9``` can be changed to ```INNER_DIAMETER = 3``` or ```INNER_DIAMETER = 25``` etc.
```
INNER_DIAMETER = 9   # Diameter of the inner "core" section
EDGE_WIDTH     = 3   # Width of the surrounding "edge" section
```
### Define Species Mixture
This program can take any number of species, however the difficulty of satisfying neighbor criteria increases with greater species number. For large species mixtures, larger core diameters will have greater chances of success.

Firstly, create a new list with the names of each species (or any taxa level) in your mixture. This list defines the mixture size and will be used for labels and colours in the exported data.
Secondly, define CHOOSE_LIST according to the name for your list.
Third, define the colours you want each species to be represented as in the exported plot image

```
LISTS = {
"SPECIES_LIST_1": ["Sorbus torminalis", "Taxus baccata", "Acer campestre", "Tilia platyphyllos"],
"SPECIES_LIST_2": ["Sorbus torminalis", "Taxus baccata", "Pinus nigra", "Prunus avium"],
"SPECIES_LIST_3": ["Acer campestre", "Tilia platyphyllos", "Pinus nigra", "Prunus avium"],
"SPECIES_LIST_4": ["Sorbus torminalis", "Taxus baccata", "Acer campestre", "Tilia platyphyllos", "Pinus nigra", "Prunus avium"]
}

CHOOSE_LIST = "SPECIES_LIST_1"

SPECIES_COLORS = {
    "Sorbus torminalis":  "#ff0000",
    "Taxus baccata":   "#020202",
    "Acer campestre":    "#ffff00",
    "Tilia platyphyllos":   "#16f138",
    "Pinus nigra":   "#0000ff",
    "Prunus avium":  "#e816f1"
}
```

### Define Criteria
MIN_NEIGHBOR_CASES decides the minimum number of cases (k) for a given species to have a certain number of neighbors (n) of a certain species. For example, if mink = 3 for n = 1 the layout must have at least 3 cases of each species being beside exactly 1 of each species. If there is any species with fewer than 3 cases where it is adjacent to exactly 1 of any species, this criterium is violated. MIN_NEIGHBOR_CASES are the more important criteria, your goal should be to maximize MIN_NEIGHBOR_CASES, especially for n that are relevent to your hypothesis. n = 0 should not be overlooked as cases where a species is no neighbors of a certain species may be used as the control treatment of a pairwise effect experiment. Be careful with high n numbers as these are progressively more difficult to satisfy and may not be relevant.

NEIGHBOR_TOLERANCE decides the allowable range of cases across species pair permutations for a given n. Rangek = maxk - mink. While this does not represent an experimental priority, limiting the range of k will discourage redundant pairwise effect representation. So, this can be used in tandem with MIN_NEIGHBOR_CASES to maximize mink for relevant n numbers. Hint: be very strict with NEIGHBOR_TOLERANCE for high n numbers. 

Criteria do not need to be met for a layout to be exported. The best layout across random seeds (described in "Define Search Effort") will be exported regardless. If the best layout meets all criteria, the following will be printed: "Exporting BEST VALID layout (all constraints met)" Otherwise, the following will be printed: "Exporting BEST OVERALL layout (constraints not fully met)"    

Below, I have given very easy criteria for a 9-diameter "core" with 4 species.

```
MIN_NEIGHBOR_CASES = {0: 1, #for n = 0, mink >= _
                      1: 1, #for n = 1, mink >= _
                      2: 1, 
                      3: 0, 
                      4: 0,
                      5: 0,
                      6: 0}

NEIGHBOR_TOLERANCE = {0: 4, #for n = 0, rangek <= _
                      1: 4, 
                      2: 4, 
                      3: 2,
                      4: 1,
                      5: 1,
                      6: 0}
```

### Define Scoring
The program works through simulated annealing optimization. First a certain number of random layouts are created, and the best one is taken. With this best one, random swaps are made between "core" trees; If the swap results in an improved layout, the swap is kept. If the swap results in a worse layout, the degree of diminishment decides the probability that the swap is kept. Why would you ever want to keep a worse swap? The program needs to be able to escape local optima, i.e. sometimes, a layout needs to get worse before it can get better. Criteria are the most important factors for comparing layouts. Scoring determines the weight that different criteria have and how the program "approaches" optimization. Different categories of points allow the incentivization of improving layouts toward criteria and improving layouts after criteria have been met. Suggestions for scoring are difficult to make. 

The following point definitions worked very well for my 9-diameter "core" plot with 4 or 6 species.

```
PERFECT_BONUS                 = 500
POINTS_PER_MIN_CONSTRAINT     = 20 # per FocalSpecies-NeighbborSpecies-minK_Success. Does not incentivize over-optimization
POINTS_PER_BALANCE_CONSTRAINT = 10 # per FocalSpecies-NeighbborSpecies-rangeK_Success. Multiplied by how much the criteria is surpassed +1

MIN_VIOLATION_PENALTY = 200 # per FocalSpecies-NeighbborSpecies-minK_Violation. Multiplied by how much the criteria is violated
TOL_VIOLATION_PENALTY = 20 # per FocalSpecies-NeighbborSpecies-rangeK_Violation. Multiplied by how much the criteria is violated.

VARIANCE_PENALTY = 2.0 # across all FocalSpecies-NeighbborSpecies case counts
RANGE_PENALTY    = 5.0 # penalizes high range case counts for a given n
ECOLOGY_WEIGHT   = 2.0 # rewards high minima of case counts (raising worst-off k for a given n), prioritzes n with low minima
MAX_EVALUATED_NEIGHBORS = 3 # Used for ecological weight ie. "raise the weakest case". Not used for MIN and TOLERANCE criteria
```
### Define Simulated Annealing Calibration
As stated in "Define Scoring", if a random swap results in a worse layout, the change in points decides the probability that the swap is kept. "Temperature" describes the probability of keeping a swap of a certain change in points (Δ). As swaps are made, temperature decreases and it becomes less likely to keep swaps that cause diminishment. Probability  ​= e^(Δ/T).
 
```
SA_T0    =  1200.0 # initial temperature, recommended: equivalent to average score difference between random swaps
SA_T_END = 0.1 #minimum temperature
```
### Define Search Effort
Search effort includes the following: 1, the number of random layouts that are generated and compared before proceeding with random swaps. 2, the number of random swaps. 3, the number of times 1 and 2 is completed, each time using a different random seed. Increasing any of these parameters will increase the number of layouts that are evaluated and compared. Keep in mind that swapping is more important for layout optimization than checking initial random layouts. For difficult criteria, I suggest you prioritize increasing the number of swaps. In addition to increasing the effort, a higher number of swaps will also decrease the rate of cooling which may improve optimization. Even for extensive searching, I recommend using a seed range of ~10 and a random attempt number no more than 10 000

Below is a very low search effort that will still successfully and quickly find a layout for the above criteria.

```
RANDOM_ATTEMPTS = 1000        # Initial random layouts
SA_SWAPS        = 5000        # Simulated annealing steps
SEED_RANGE      = range(2)    # Number of random seeds
```
### Define Export
The program will create four export files. A few input parameters will automatically be included in the exported file names. You can decide if "Test" is added to the start of these file names. "Test" may be useful when exploring parameters as you will be creating many files that you may later want to delete. You must define the global location of the folder where you want the files to be saved. If a file with the same name already exists, a numeric suffix is added automatically.

1, a vector image of the plot in .svg format. 

2, an Excel spreadsheet with the counts of cases for each n for each species-pair permutation in .xlsx format. It also includes summary statistics with the input criteria parameters and the minimums and ranges in the exported layout. These exported mins and ranges are colour coded according to if their attempted criteria were failed -> yellow, met -> green, or surpassed -> blue. 

3, an Excel spreadsheet with a summary of the exported layout. It records each tree's position, species, and neighbor diversity metrics. Position metrics include hexagonal coordinates, "zone" defined as whether the tree is in the core or edge section, and "hex radius" defined as the number of hex-rings from the plot center. Neighbor diversity metrics evaluate the composition of that tree's immediate neighbors. Diversity metrics include species richness, Shannon index, Pielou's evenness index, and the identity of the most represent neighbor (multiple are written in the case of a tie). Row order is spirally organized with the center tree at the top of the spreadsheet, followed by the rightmost tree in the first hex-ring, followed by the clockwise sequence of trees in the same ring, followed by the rightmost tree in the next hex-ring, etc.

4, a text document with each input parameter in .txt format. This can be useful for parameter exploration and reproducibility.  

```
TEST = True # If True, adds 'Test' prefix to output filenames.
EXPORT_DIR = r"C:\MyExperiment\Method\PlantingLayout\AlphaPlot\Exports\c9e3\Tests" # Where do you want the exports to be saved? 
```

### Run the script

From the command line:
```
python AlphaPlot.py
```
Or run directly from an IDE (e.g. VS Code, PyCharm).

## Future Improvements

For suggestions or bug reports, feel free to contact me via email: asc.sheppard@gmail.com 

## AI Disclosure

I received significant help from ChatGPT for coding and designing this program. The above workflow was written by myself with minimal input from AI.

## Reference List

Jerome K Vanclay, Experiment designs to evaluate inter- and intra-specific interactions in mixed plantings of forest trees, Forest Ecology and Management, Volume 233, Issues 2–3, 2006, Pages 366-374, ISSN 0378-1127, https://doi.org/10.1016/j.foreco.2006.05.034

Bruelheide, H., Nadrowski, K., Assmann, T., Bauhus, J., Both, S., Buscot, F., Chen, X., Ding, B., Durka, W., Erfmeier, A., Gutknecht, J. L. M., Guo, D., Guo, L., Härdtle, W., He, J., Klein, A., Kühn, P., Liang, Y., Liu, X., Michalski, S., Niklaus, P., Pei, ., Scherer-Lorenzen, M., Scholten, T., Schuldt, A., Seidler, G., Trogisch, S., von Oheimb, G., Welk, E., Wirth, C., Wubert, T., Yang, X., Yu, M., Zhang, S., Zhou, H., Fischer, M., Ma, K., & Schmid, B. (2014). Designing forest biodiversity experiments: general considerations illustrated by a new large experiment in subtropical China. Methods in Ecology and Evolution, 5(1), 74–89. https://doi.org/10.1111/2041-210X.12126
