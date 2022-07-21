library(torch)


GraphConvolution <- nn_module(
  initialize = function(in_features, out_features, bias=TRUE){
    # torch_manual_seed(123)
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
    print("self$stdv is ")
    print(self$stdv)
    print("self$bias is")
    print(self$bias)
    print("self$weight is")
    print(self$weight)
    # if (!is.na(self$bias)) {
      self$bias <- with_no_grad(self$bias$data()$uniform_(-stdv, stdv))
    # }
  },
  forward = function(input, adj){
    support <- torch_mm(input, self$weight)
    output <- torch_mm(adj, support)
    # if (!is.na(self$bias)){
      print("output+self$bias")
      print(output+self$bias)
      return(output+self$bias)
    # }else {
      # return(output)
    # }
  }
)

# a = torch_tensor(matrix(c(1,2,3,4,5,6),nrow = 3))
# b = torch_tensor(matrix(c(1,2,3,4,5,6,7,8,9),nrow = 3))
# gc = GraphConvolution(2,2)
# dec = gc(a,b)








