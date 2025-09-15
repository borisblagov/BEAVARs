Current main dev branch is "testingDispatch", which is the direction I want to go

To do
- [ ] Use multiple dispatch constructors to allocate output after test dispatch for forecasting and IRFs


## 15.09.2025

Bug: You cannot have NaNs in the low-frequency file, this would throw an error. E.g. if your quarterly data looks like this

|data| Y |
|--- |---|
| Q1 | 23.0|
| Q2 | NaN|

you will get an error. This has to be fixed so that one can use low-frequency data with a ragged edge
