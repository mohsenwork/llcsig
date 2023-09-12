# main.R

# Import the directed_reachable function from utilities.R
# Get the current working directory
# cwd <- getwd()

# Print the current working directory
source("./simulations/common/R/utilities.R")
library(reticulate)


main <- function() {
  # Get the command-line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check if the required number of arguments is provided
  if (length(args) < 1) {
    stop("Insufficient number of arguments.")
  }
  
  data_path = args[1]
#   print(data_path)
  df = py_load_object(data_path, pickle = 'pickle')

    # print("R Data:")
    # for (name in names(df)) {
    #     print(name)
    #     print(df[[name]])
    # }
    # print("R End Data")

  # Call the directed_reachable function with the provided arguments
  result <- directed_reachable(df[['x']], df[['y']], df[['C']], unlist(df[['J']]), df[['D']], df[['B']], df[['U']])
  
  print(result)
}

# Call the main function
main()
