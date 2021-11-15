#' Variable Selection
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#' @importFrom purrr quietly
#'
#' @examples
var_select_cat <- function(x,y){
  # 自动会删去那些为0的level
  if(anyNA(x)){
    levels(x) = c(levels(x), 'MissingDalian')
    x[is.na(x)] <- 'MissingDalian'
  }
  fit = purrr::quietly(chisq.test)(x,y)$result
  # ans = ifelse(fit$parameter > 1L, wilson_hilferty(fit$statistic, fit$parameter), fit$statistic)
  # 在不小于2.2E-16的时候，先用原本的p-value。
  ans = ifelse(fit$parameter > 1L, ifelse(fit$p.value > 10^(-16),
                                          qchisq(1-fit$p.value, df = 1),
                                          wilson_hilferty(fit$statistic,fit$parameter)),fit$statistic)
  return(ans)
}

wilson_hilferty = function(chi,df){ # 把df = K 的卡方变成 df = 1 的卡方
  ans = max(0, (7/9 + sqrt(df) * ( (chi / df) ^ (1/3) - 1 + 2 / (9 * df) ))^3)
  return(ans)
}

# 目前改成分位数法
var_select_noncat <- function(x, y, Nt, Jt){
  # Followed Loh09 Paper
  m = mean(x,na.rm = T)
  s = sd(x,na.rm = T)
  if(Nt >= 30 * Jt){
    split_1 = quantile(x,c(0.25,0.5,0.75), na.rm = T)
    if(length(unique(split_1)) != 3){
      split_1 = c(m - s *sqrt(3)/2, m, m + s *sqrt(3)/2)
    }
    x = cut(x, breaks = c(-Inf, split_1, Inf), right = TRUE)
    # 这里如果 right = FALSE, 将会违背我们分位数的初衷
  }else{
    split_2 = quantile(x,c(0.33,0.66), na.rm = T)
    if(length(unique(split_2)) != 2){
      split_2 = c(m - s *sqrt(3)/3, m + s *sqrt(3)/3)
    }
    x = cut(x, breaks = c(-Inf, split_2, Inf), right = TRUE)
  }
  return(var_select_cat(x,y))
}

# 被淘汰的方法
# var_select_noncat <- function(x, y, Nt, Jt){
#   # Followed Loh09 Paper
#   m = mean(x)
#   s = sd(x)
#   if(Nt >= 30 * Jt){
#     x = cut(x, breaks = c(-Inf, m - s *sqrt(3)/2, m, m + s *sqrt(3)/2, Inf), right = FALSE)
#   }else{
#     x = cut(x, breaks = c(-Inf, m - s *sqrt(3)/3, m + s *sqrt(3)/3, Inf), right = FALSE)
#   }
#   return(var_select_cat(x,y))
# }


var_select_all <- function(x, y, Nt, Jt){
  if(class(x) %in% c('numeric', 'integer')){
    return(var_select_noncat(x, y, Nt, Jt))
  }else{
    return(var_select_cat(x,y))
  }
}

# 这个函数是为了防止LDA因为constant in groups 而报错
within_check <- function(y,x){
  # 进来的数据应该没有NA
  # 返回TRUE，就代表着这个数据within group constant
  if(class(x) %in% c('numeric', 'integer')){
    x = round(x,8) # 有些e-16的变化都被R察觉到了，但是LDA察觉不到
    return(within_check_helper(y,x))
  }else{
    if(length(unique(x) == 1)){
      return(TRUE)
    }
    mmm = model.matrix(y~x)[,-1,drop = FALSE] # LDA好像不关心第一列
    res = apply(mmm,2,function(x_x) within_check_helper(y,x_x))
    return(any(res))
  }
}


within_check_helper <- function(y,x){
  overall_list = unique(y)
  for(i in 1:length(overall_list)){ #对于每一个y下面的subgroup
    idx_tmp = which(y == overall_list[i])
    if(length(unique(x[idx_tmp])) > 1){ # 如果X只要有一个非constant，就能估计方差
      return(FALSE)
    }
  }
  return(TRUE)
}

naive_impute <- function(x){
  for(i in 1:ncol(x)){
    if(anyNA(x[,i])){
      if(class(x[,i]) %in% c('numeric', 'integer')){
        x[which(is.na(x[,i])),i] = mean(x[,i], na.rm = T)
      }else{
        x[which(is.na(x[,i])),i] = getmode(x[,i])
      }
    }
  }
  return(x)
}

getmode <- function(v) {
  uniqv <- unique(v)
  return(uniqv[which.max(tabulate(match(v, uniqv)))])
}

best_friend <- function(x_new, x_original){ # 找到最像的朋友们，按相似度从高到底排序
  # check if they have the same dimension
  # 这里以后补充一个get_error function for robustness
  # 还要加一个new level的处理方法

  idx_keep = which(!is.na(x_new))
  x_new_tmp = x_new[idx_keep]
  x_original_tmp = x_original[idx_keep]
  # 删去原数据中NA的那些行
  idx_notNA = which(apply(x_original_tmp,1, function(x) !any(is.na(x))))
  x_original_tmp = x_original_tmp[idx_notNA,]
  # data check
  # 马氏距离：对于factor，变成了one-hot-encoding之后，cov矩阵不可逆
  # 以后采用标准化之后的欧氏距离
  # 为了防止dummyX使得biased towards factor，我们先变一下, binary
  x_combined = rbind(x_new_tmp,x_original_tmp)
  for(i in 1:ncol(x_combined)){
    if(!class(x_combined[,i]) %in% c('numeric', 'integer')){
      x_combined[,i] = (x_combined[,i] == x_combined[1,i]) + 0
    }
  }
  x_combined = scale(x_combined,center = TRUE,scale = TRUE)
  dist_tmp = apply(x_combined,1,function(x) dist(rbind(x_combined[1,],x)))
  return(idx_notNA[order(dist_tmp[-1])]) # 返回一个距离的排序
}









