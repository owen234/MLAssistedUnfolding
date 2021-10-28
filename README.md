# MLAssistedUnfolding

Input files can be found in this directory at DESY: ```/nfs/dust/h1/group/owen/MLAssistedUnfolding-files ```

The Jupyter notebooks run the NN regression training (**NN_h1_reg_v4b.ipynb**) and making an ntuple from the regression DNN (**NN_h1_make_dnn_tree_v2.ipynb**)

The code **fill_hists1** reads in the DNN output TTree and makes input plots for the unfolding.

To run it, do something like this
```
<start root>
.L fill_hists1.c
fill_hists f("dnn-output-file-name.root")
f.Loop()
```

The unfolding code (**unfold1.c** and **unfold_comp1.c**) used tutorials/unfold/testUnfold1.C from the 6.24 ROOT distribution as a template.

