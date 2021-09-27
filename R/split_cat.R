#' Splitting in the LDAtree
#'
#' @param x
#' @param y
#' @param mis_curr
#'
#' @return
#' @export
#' @importFrom MASS lda
#'
#' @examples
split_cat <- function(x,y, mis_curr){
  x = factor(x) # 这个relevel很有必要，要不然会因为model.matrix有一列全是0从而导致constant
  fit = MASS::lda(y~x-1)
  l = levels(x)[which(coefficients(fit)>0)] # 无法判断是否有提升
  return(l)
}

split_noncat <- function(x,y,datx, mis_curr){ # 这一步跑得太太太太慢
  threshold = sort(unique(x))
  ans = rep(Inf,length(threshold))
  for(i in 1:(length(threshold)-1)){
    idx = which(x<= threshold[i])
    y_l = y[idx]
    y_r = y[-idx]
    mis_l = mis_r = 0
    if(length(unique(y_l))!=1){ # 我们采用小于等于的准则
      x_l = datx[idx,]
      # fit = quietly(total_LDA)(x_l, y_l)$result # 用子节点的LDA结果作为划分的依据
      # mis_l = sum(predict(fit,cbind(x_l,y_l))$class != y_l)
      mis_l = get_error_LDA(x_l, y_l)
    }
    if(length(unique(y_r))!=1){ # 我们采用小于等于的准则
      x_r = datx[-idx,]
      # fit = quietly(total_LDA)(x_r, y_r)$result
      # mis_r = sum(predict(fit,cbind(x_r,y_r))$class != y_r)
      mis_r = get_error_LDA(x_r, y_r)
    }
    ans[i] = mis_l + mis_r
  }
  idx_threshold = which(ans == min(ans))[1]
  # print(c(ans[idx_threshold],mis_curr))
  if(ans[idx_threshold] >= mis_curr){
    return(NULL)
  }else{
    return(threshold[idx_threshold])
  }
}
