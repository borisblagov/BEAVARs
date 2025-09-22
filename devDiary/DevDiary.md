Current main dev branch is "testingDispatch", which is the direction I want to go

To do
- [ ] Use multiple dispatch constructors to allocate output after test dispatch for forecasting and IRFs


## 15.09.2025

OPEN Bug: You cannot have NaNs in the low-frequency file, this would throw an error. E.g. if your quarterly data looks like this

|data| Y |
|--- |---|
| Q1 | 23.0|
| Q2 | NaN|

you will get an error. This has to be fixed so that one can use low-frequency data with a ragged edge


## 22.09.2025

- [ ] in branch testingDispatchNewNew add a new function for the intercept location. Either like makeHypSetup(), where it selects it based on the model type or like selectModel() where it takes the string
- [ ] Currently the functions dispatchModel() are in their respective files. The main beavar() function for CPZ and Chann2020minn is still in the main file. Decide where you want them
- [ ] in Chan2020minn think whether to add new makeDataSetup for non Time-Arrays