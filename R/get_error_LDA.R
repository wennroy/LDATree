#' LDA Related Functions
#'
#' Calculate the error generated from LDA when fitting the nodes.
#'
#' @param x
#' @param y
#'
#' @return
#'
#' @importFrom MASS lda
#' @export
#'
#' @examples
get_error_LDA <- function(x,y,prior){
  # 两种 LDA 打包在一起
  # 函数的目的是输入x和y，换回一个错误的个数
  # 分类方式有LDA，众数以及constant group的特殊方法
  # flag_lda = TRUE
  if(is.null(dim(x))){
    x = matrix(x,,1)
  }
  x = naive_impute(x) # 把NA先填上，到这一步的时候不会出现一列全是NA的情况
  x = droplevels(x) # 不用的level去掉
  dat_lda = as.data.frame(cbind(x,y))

  mis_class = tryCatch({
    fit <<- MASS::lda(y~., data = dat_lda, prior = prior) # fit an LDA model in the node
    return(sum(predict(fit,dat_lda)$class != y))
  }, error = function(e) { # variable being constant within groups，应该只有这一种错误
    # flag_lda <<- FALSE

    ### 考虑prior进去，来算众数
    tab_tmp = tabulate(y) * prior
    mis_class_mode = sum(tabulate(y)[-which.max(tab_tmp)])
    ###

    ### 开始复杂的函数
    idx_constant = which(apply(x,2,function(x) within_check(y,x))) # 那些constant的组
    for(i in idx_constant){
      if(class(x[,i]) %in% c('numeric', 'integer')){
        x[,i] = x[,i]+rnorm(dim(x)[1],sd = mean(x[,i])/100) # 对于numeric的变量，加个小数使得能求方差

      }
    }



    return() # 众数估计
  })

  # 开始处理constant group
  x = x[,, drop = FALSE] # 改掉group constant
  # x = x[,apply(x,2,function(t) length(unique(t))!=1)] # 讲那些列一样的去掉，这一步主要是在防splitting
  return(mis_class)
}


# 这个函数的目的是对叶子结点进行预测
# 并返回一个可以用来预测其他变量的function
# pred_LDA <- function(x,y,prior){
#   flag_lda = TRUE
#   m = 'lda'
#
#   ### 把x变成一个矩阵
#   if(is.null(dim(x))){
#     x = matrix(x,,1)
#   }
#   ###
#
#   x = naive_impute(x) # impute NAs
#   # x = x[,apply(x,2,function(x) !within_check(y,x)), drop = FALSE] # 改掉group constant
#   # x = droplevels(x)
#   # x = x[,apply(x,2,function(t) length(unique(t))!=1)] # 讲那些列一样的去掉，这一步主要是在防splitting
#   result = tryCatch({
#     fit <<- MASS::lda(y~., data = x, prior = prior) # fit an LDA model in the node
#   }, error = function(e) { # variable being constant within groups
#     flag_lda <<- FALSE
#   })
#   if(!flag_lda){ # 用众数来估计
#     m = 'mode'
#     tab_tmp = tabulate(match(y, unique(y)))
#     ans = unique(y)[which.max(tab_tmp)]
#   }else{
#     ans = fit
#   }
#   return(list(m,ans))
# }



pred_LDA <- function(x,y,prior){

  ### 把x变成一个矩阵
  if(is.null(dim(x))){
    x = matrix(x,,1)
  }
  ###

  ### 考虑prior进去，来算众数
  # m = 'mode'
  tab_tmp = tabulate(y) * prior
  ans_mode = levels(y)[which.max(tab_tmp)]
  mis_class_mode = sum(y != ans_mode) # misclass = 1的简单情况
  ###

  x = naive_impute(x) # impute NAs
  # x = x[,apply(x,2,function(x) !within_check(y,x)), drop = FALSE] # 改掉group constant
  # x = droplevels(x)
  result = tryCatch({
    fit <<- MASS::lda(y~., data = dat_lda, prior = prior) # fit an LDA model in the node
    ans = fit
    return(list('lda',ans,sum(predict(fit,dat_lda)$class != y)))
  }, error = function(e) { # variable being constant within groups，应该只有这一种错误
#
#     ### 开始复杂的函数
#     idx_constant = which(apply(x,2,function(o_o) within_check(y,o_o))) # 那些constant的组
#     idx_cat = c()
#     for(i in idx_constant){
#       if(class(x[,i]) %in% c('numeric', 'integer')){
#         x[,i] = x[,i]+rnorm(dim(x)[1],sd = mean(x[,i])/100) # 对于numeric的变量，加个小数使得能求方差
#       }else{
#         idx_cat = c(idx_cat,i) # 将categorical var 存起来
#       }
#     }
#     if(length(idx_cat) == 0){ # 如果全是numeric variable
#       fit <<- MASS::lda(y~., data = dat_lda, prior = prior)
#       ans = fit
#       return(list('lda',ans))
#     }else{
#       # 产生一个分类的函数
#       table_tmp = unique(cbind(x[,idx_cat],y)) # 记录了每一个X唯一对应的Y是什么
#       Jt = length(unique(y))
#       if(dim(unique(table_tmp[,1:length(idx_cat)]))[1] == Jt){ # 如果能全部覆盖的话，就直接100%预测了
#         # 进行一个直接one-to-one的mapping
#       }else{
#         # 进行一个one-to-one的mapping + 其他不能预测的再跑一个LDA
#         # 先看有哪些y被唯一的表达了
#         record_x = apply(table_tmp[,-length(idx_cat)-1],1, function(o_o) paste(o_o,collapse = ''))
#         row_unique = which(record_x %in% names(table(record_x))[which(table(record_x) == 1)])
#         # 先进行一个one-to-one mapping
#         result_f <- function(x){
#         }
#       }
#       # table_unique = table_tmp[table_tmp[,1] %in% which(table(table_tmp[,1])==1),]
#       # # 记录下来那些可以准确预测Y的X们
#       # record_list = rbind(record_list,cbind(i,table_unique))
#     }
#
#     return() # 众数估计
    return(list('mode', ans_mode, mis_class_mode))
  })
  # return(ifelse(result[[2]] < ans_mode, result, list('mode', ans_mode)))
  return(result)
}
