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
  m = mean(x)
  s = sd(x)
  if(Nt >= 30 * Jt){
    split_1 = quantile(x,c(0.25,0.5,0.75))
    if(length(unique(split_1)) != 3){
      split_1 = c(m - s *sqrt(3)/2, m, m + s *sqrt(3)/2)
    }
    x = cut(x, breaks = c(-Inf, split_1, Inf), right = TRUE)
    # 这里如果 right = FALSE, 将会违背我们分位数的初衷
  }else{
    split_2 = quantile(x,c(0.33,0.66))
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
