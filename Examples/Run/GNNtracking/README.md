Implement L2IT GNN development into ACTS
========================================
Repository summarizing the different steps for a novel technique to a full track reconstruction using graph neural networks


Table of contents
=================

  * [1 Compile ACTS](#1-compile-acts)
    * [1.1 Build and compile](#11-build-and-compile)
  * [2 Generate and simulate events](#2-generate-and-simulate-events)
  * [3 Module map creation](#3-module-map-creation)



# 1 Compile ACTS

# 1.1 Build and compile

Let's get started!

```
cd acts
mkdir build
source CI/setup_cvmfs_lcg.sh
cmake -B build -S . -DACTS_BUILD_FATRAS=on -DACTS_BUILD_EXAMPLES=on -DACTS_BUILD_EXAMPLES_PYTHIA8=on -DACTS_BUILD_PLUGIN_TGEO=on  -DACTS_BUILD_GNN_TRACKING=on
cmake --build build -j 4 ## this number can be changed, it represents the number of threads one would compile on
```

# 2 Module map creation

To create the module map:

```

build/bin/ActsRecGNNTracks \
    --input-dir=INPUT_DIR_DATA \
    --output-dir=OUTPUT_DIR_MODULE_MAP \
    --min-pt-cut=0.5 \
    --min-nhits=3 \
    --root-filename=ROOT_FILE_NAME.root \
    --give-cut-values=True \

```

With:
  * input-dir: the repository containing the data (only in .csv format for now)
  * output-dir: the repository to save the module map.
  * min-pt-cut: the pt cut applied on simulated particles
  * min-nhits: the minimum number of hits a simulated particles must have
  * root-filename: the name of the root file containing the module map
  * give-cut-values: boolean to add cut values on {z0,phislope,deta,dephi} for each module connections

# 3 Graph Creation

  To create graphs:

```

build/bin/ActsRecGNNTracks \
    --input-dir=/sps/l2it/CommonData/event_data/ttbar_PU200_tuningATLAS/1event/ \
    --input-module-map=/sps/l2it/crougier/GitLab/l2it_acts/GraphCreation/GraphCreationTools_1000event.root \
    --output-dir=/sps/l2it/crougier/GitLab/l2it_acts/GraphCreation/graph/ \
    --give-true-graph=False \
    --save-graph-on-disk=False \
    --min-pt-cut=0.5 \
    --min-nhits=3 \

```

With:
  * input-dir: repository of input event data
  * input-module-map: input module map
  * output-dir: where to save the graph if needed
  * give-true-graph: boolean, create true graphs on request. Default to false.
  * save-graph-on-disk: boolean, save graphs on disk. Default to false. 
  * min-pt-cut: the pt cut applied on simulated particles
  * min-nhits: the minimum number of hits a simulated particles must have


    All configurations given in example are the one used to produce all results presented in: https://arxiv.org/abs/2103.00916


