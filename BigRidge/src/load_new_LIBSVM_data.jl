
include("dataLoad.jl")
initDetails()

datasets = ["madelon.t"] #  w1a, SUSY,
for  dataset in datasets
transformDataJLD(dataset)
X,y = loadDataset(dataset) #
showDetails(dataset)
end
