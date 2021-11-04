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

I have tried **RooUnfold** in addition to **TUnfold**.  I found that RooUnfoldBayes gives nonsense output if I use the default initial prior, which is to set it to the true distribution.  It basically converges quickly to the true distribution and has tiny unfolding uncertainties for the bins with the worst resolution and the correlation matrix for those bins is all close to 1.  If I instead initialize it to a flat prior and run it for MANY iterations (several thousand), it converges to the same output as TUnfold.  The hack for RooUnfoldBayes to start with a flat prior is in the few lines below.

```
void RooUnfoldBayes::setup()
{
...
  if (_N0C!=0.0) {
    //*** Owen: try a flat prior to start.
    //_P0C= _nCi;
    //_P0C *= 1.0/_N0C;
    _P0C = 1.0/_nc ;
  }
}
```
