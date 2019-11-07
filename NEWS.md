# DiversityOccupancy 1.1.0

* Added a `NEWS.md` file to track changes to the package.
* `batchoccu()` now skips species that don't converge and shows which species are those. 
* `batchoccu()` gains a `SppNames` parameter to name the species models and graphs.
* `diversityoccu()` is no longer dependant on the `glmmulti` package, and now the model selection is dependant on `MuMIn`
* `responseplot.occu()` Now can handle categorical variables
