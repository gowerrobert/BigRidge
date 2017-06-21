
include("dataLoad.jl")
initDetails()

datasets = ["protein"] #  w1a, SUSY,
for  dataset in datasets
transformDataJLD(dataset)
X,y = loadDataset(dataset) #
showDetails(dataset)
end
