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
get_error_LDA <- function(x,y){
  # 两种 LDA 打包在一起
  flag_lda = TRUE
  if(is.null(dim(x))){
    x = matrix(x,,1)
  }
  x = x[,apply(x,2,function(t) length(unique(t))!=1)] # 讲那些列一样的去掉，这一步主要是在防splitting
  dat_lda = as.data.frame(cbind(x,y))
  result = tryCatch({
    fit <<- MASS::lda(y~., data = dat_lda) # fit an LDA model in the node
    mis_class = sum(predict(fit,dat_lda)$class != y)
  }, error = function(e) { # variable being constant within groups
    flag_lda <<- FALSE
  })
  if(!flag_lda){ # 用众数来估计
    # 这个quietly好强，直接干掉了所有的东西
    # s = quietly(stepclass)(y~., data = dat_lda, method = 'lda', direction = 'forward', fold = 10, criterion = "AS")$result # 目前这个函数十分不稳定，我觉得日后得抛弃
    # fit <- lda(s$formula, data = dat_lda) # 把得到的formula再放进lda跑一遍
    tab_tmp = tabulate(match(y, unique(y)))
    mis_class = sum(tab_tmp) - max(tab_tmp) # 众数估计
  }
  return(mis_class)
}
# stepclass(x_l, y, method = 'lda', direction = 'forward', fold = 10)

# 把前面的函数改一改
pred_LDA <- function(x,y){
  flag_lda = TRUE
  m = 'lda'
  if(is.null(dim(x))){
    x = matrix(x,,1)
  }
  x = x[,apply(x,2,function(t) length(unique(t))!=1)] # 讲那些列一样的去掉，这一步主要是在防splitting
  result = tryCatch({
    fit <<- MASS::lda(y~., data = x) # fit an LDA model in the node
  }, error = function(e) { # variable being constant within groups
    flag_lda <<- FALSE
  })
  if(!flag_lda){ # 用众数来估计
    m = 'mode'
    tab_tmp = tabulate(match(y, unique(y)))
    ans = unique(y)[which.max(tab_tmp)]
  }else{
    ans = fit
  }
  return(list(m,ans))
}
