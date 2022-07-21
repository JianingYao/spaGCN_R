library(torch)


GraphConvolution <- nn_module(
  initialize = function(in_features, out_features, bias=TRUE){
    self$in_features = in_features
    self$out_features = out_features
    self$weight = nn_parameter(torch_empty(in_features, out_features))
    self$bias = nn_parameter(torch_empty(out_features))
    self$reset_parameters()
  },
  reset_parameters = function(){
    stdv <- 1/sqrt(self$weight$size(2))
    self$stdv <- stdv
    self$weight <- with_no_grad(self$weight$data()$uniform_(-stdv, stdv))
    self$bias <- with_no_grad(self$bias$data()$uniform_(-stdv, stdv))
  },
  forward = function(input, adj){
    support <- torch_mm(input, self$weight)
    output <- torch_mm(adj, support)
    return(output+self$bias)
  }
)









