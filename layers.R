library(torch)


reset_parameters <- function(gc){
  stdv <- 1/sqrt(gc$weight$size(1))
  gc$stdv <- stdv
  gc$weight <- with_no_grad(gc$weight$uniform_(-stdv, stdv))
  if (!is.na(gc$bias)) {
    gc$bias <- with_no_grad(gc$bias$uniform_(-stdv, stdv))
  }
}


GraphConvolution <- nn_module(
  initialize = function(in_features, out_features, bias=TRUE){
    self$in_features = in_features
    self$out_features = out_features
    self$weight = nn_parameter(torch_empty(in_features, out_features))
    self$bias = nn_parameter(torch_empty(out_features))
    reset_parameters(self)
  },
  forward = function(input, adj){
    support <- torch_mm(input, self$weight)
    output <- torch_mm(adj, support)
    if (!is.na(self$bias)){
      return(output+self$bias)
    }else {
      return(output)
    }
  }
)



# 
# GraphConvolution <- nn_module(
#   reset_parameters.GraphConvolution = function(){
#     stdv <- 1/sqrt(self$weight$size(1))
#     self$weight <- with_no_grad(self$weight$uniform_(-stdv, stdv))
#     if (!is.na(self$bias)) {
#       self$bias <- with_no_grad(self$bias$uniform_(-stdv, stdv))
#     }
#   },
#   initialize = function(in_features, out_features, bias=TRUE){
#     self$in_features = in_features
#     self$out_features = out_features
#     self$weight = nn_parameter(torch_empty(in_features, out_features))
#     self$bias = nn_parameter(torch_empty(out_features))
#     reset_parameters()
#   },
#   forward = function(input, adj){
#     support <- torch_mm(input, self$weight)
#     output <- torch_mm(adj, support)
#     if (!is.na(self$bias)){
#       return(output+self$bias)
#     }else {
#       return(output)
#     }
#   }
# )











