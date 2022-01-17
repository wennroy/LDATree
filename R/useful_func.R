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


var_select_all <- function(x, y, Nt, Jt, select_method){
  # 把方法外包出去
  if(pmatch(select_method, c('LDATree', 'FACT')) == 1){
    # LDATree
    return(var_select_LDATree(x, y, Nt, Jt))
  }else if(pmatch(select_method, c('LDATree', 'FACT')) == 2){
    # FACT
    return(var_select_FACT(x, y, Nt, Jt))
  }
}

var_select_LDATree <- function(x, y, Nt, Jt){
  if(class(x) %in% c('numeric', 'integer')){
    return(var_select_noncat(x, y, Nt, Jt))
  }else{
    return(var_select_cat(x,y))
  }
}

var_select_FACT <- function(x, y, Nt, Jt){
  # if(class(x) %in% c('numeric', 'integer')){
  #   F_stat = anova(lm(x~y))$`F value`[1] # F-value
  #   return(F_stat)
  # }else{
  #   return(var_select_cat(x,y))
  # }
  F_stat = anova(lm(x~y))$`F value`[1] # F-value
  return(F_stat)
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

# FACT cat -> num

# x <- dat[,o_o]
# y = response

fact_cat <- function(x,y, prior){
  dummy_matrix = model.matrix(y~x-1) # Get the dummy matrix
  # fit = eigen(cov(dummy_matrix)) # Eigen decomposition, 为了防止LDA矩阵不可逆
  # eigen_keep = which(round(fit$values,8) > 0) # 保留正值
  # X_dummy = dummy_matrix %*% fit$vectors[,eigen_keep] # Projection
  fit = princomp(dummy_matrix, cor = TRUE)
  X_dummy = fit$scores[,round(fit$sdev^2,8) > 0, drop = FALSE] # 这个drop = F感觉得变成default，总会报错
  X_dummy = X_dummy[,apply(X_dummy,2,function(x_x) !within_check(y,x_x)), drop = FALSE] # 改掉group constant
  new_data = data.frame(y, X_dummy)
  fit_lda = lda(y~., data = new_data, prior = prior) # 这个LDA需要后面的scaling，所以先不变robust了
  X_num = X_dummy %*% fit_lda$scaling[,1] # Project 到 LD1 上面去
  reexp = unique(data.frame(x,X_num)) # 找到x和X_num的对照表
  return(list(X_num,reexp))
}

PCA_Y <- function(datx, response_tmp, beta_ratio = 0.05, situation = 1){
  # 参数这么多，完全是为了prediction
  d_matrix = as.matrix(datx)
  fit = princomp(d_matrix, cor = TRUE)
  # 这里需要记录下来转移的均值和系数，供之后的prediction使用
  comp_kept = (fit$sdev^2 >= (fit$sdev^2)[1] * beta_ratio)
  X_d = fit$scores[,comp_kept] # 这个决定了dat_tmp_y有多少列

  # 全在这里面再算一遍，robustness，因为居然出现了y和w维数不同的错误
  # 确实不明白是怎么发生的，但是全在这里面算一遍有助于解决问题
  trans_w_tmp = scale(X_d,center = TRUE, scale = FALSE)
  dat_tmp_w = abs(trans_w_tmp)
  p_stat = p.adjust(apply(dat_tmp_w,2,function(o_o) car::leveneTest(y = o_o, group = response_tmp)[[3]][1]))
  levene_sig <- (p_stat <= 0.05)

  # Function to get y from x
  linear_split_trans <- function(x_new){
    y_new = matrix((as.matrix(x_new) - fit$center) / fit$scale,1) %*% as.matrix(fit$loadings[,comp_kept])
    return(y_new)
  }
  if(situation == 3){
    linear_split_trans <- function(x_new){
      y_new = matrix((as.matrix(x_new) - fit$center) / fit$scale,1) %*% as.matrix(fit$loadings[,comp_kept])
      w_new = abs(y_new - attr(trans_w_tmp,'scaled:center'))
      return(w_new)
    }
    return(linear_split_trans)
  }else if(situation == 4){
    linear_split_trans <- function(x_new){
      y_new = matrix((as.matrix(x_new) - fit$center) / fit$scale,1) %*% as.matrix(fit$loadings[,comp_kept])
      w_new = abs(y_new - attr(trans_w_tmp,'scaled:center'))
      w_new = nsphere(w_new)
      return(w_new)
    }
    return(linear_split_trans)
  }else if(situation == 5){
    linear_split_trans <- function(x_new){
      y_new = matrix((as.matrix(x_new) - fit$center) / fit$scale,1) %*% as.matrix(fit$loadings[,comp_kept])
      w_new = abs(y_new - attr(trans_w_tmp,'scaled:center'))
      w_new[,levene_sig] = nsphere(w_new[,levene_sig])
      return(w_new)
    }
    return(linear_split_trans)
  }
  return(list(X_d, linear_split_trans))
}

############
# 浓缩成一个函数
# 输出结果包括：是否继续划分，如果是的话，返回cut point, idx_children, no_class

fact_univar <- function(x, y, node_tmp, prior, misclass_cost, Nj, min_nsize, get_size, simple_mean = FALSE){
  if(!simple_mean){
    # prior calculation
    pjt = prior * node_tmp$portion / Nj
    pjgt = pjt / sum(pjt) # standardize
    # node_error_rate = sum(pjt) - max(pjt) # 这里要用绝对错误，而不是相对错误
    node_error_rate = mis_cost_cal(proportion = pjt, misclass_cost = misclass_cost)
    threshold = split_fact_uni(x,y, pjgt, misclass_cost, min_nsize)
    if(is.null(threshold)){
      # node_saved[[node_tmp$idx]] = list(node_tmp)
      return(NULL) # 停止分割
    }
    # 下面可以继续分
    node_tmp$split_cri = round(threshold,1) # 只保留一位小数
  }else{
    node_error_rate = Inf
    node_tmp$split_cri = mean(x) # Take mean to prevent early stopping
  }

  no_class = length(node_tmp$split_cri) + 1 #  产生了几个子节点
  group_idx = cut(x,breaks = c(-Inf,node_tmp$split_cri,Inf),
                  labels = 1:no_class,right = TRUE) # 包含右边界
  idx_children = sapply(1:no_class, FUN = function(o_o) which(group_idx == o_o)) # save the idx_r for children
  # subnode_index_c = node_tmp$idx_c # 剩下的cov

  # 出口4.5: 判断是否不存在任何一组大于minimum node size，到达了便退出
  # print(no_class)
  # print(node_tmp$split_cri)
  ephemeral = sapply(seq(no_class), function(o_o) sum(table(y[idx_children[[o_o]]]) >= min_nsize))
  # print(ephemeral)
  if(min(ephemeral) < 1){
    # print(table(response_tmp[idx_children[[o_o]]]))
    # cat('NA happens there 1:', node_tmp$idx,'\n')
    print('Exit 4.5')
    # node_tmp$split_idx = node_tmp$split_cri = NA
    # node_saved[[node_tmp$idx]] = list(node_tmp)
    return(NULL) # 停止分割
  }


  # 出口5: 决定是否要继续分了，根据错误率
  # 只有pre-stopping在乎
  if(pmatch(get_size, c('CV', 'pre-stopping')) == 2){
    children_error_rate = numeric(no_class)
    for(o_o in 1:no_class){
      pjt = prior * table(y[idx_children[[o_o]]]) / Nj
      # pjgt = pjt / sum(pjt)
      # children_error_rate[o_o] = sum(pjt) - max(pjt)
      children_error_rate[o_o] = mis_cost_cal(proportion = pjt, misclass_cost = misclass_cost)
    }

    if(sum(children_error_rate) >= node_error_rate){
      print('Exit 5')
      # node_tmp$split_idx = node_tmp$split_cri = NA
      # node_saved[[node_tmp$idx]] = list(node_tmp)
      return(NULL) # 停止分割
      # next # 出口5
    }
  }
  return(list(node_tmp, idx_children, no_class))
}

############

# 这个函数的目的是给定各个组的proportion
# 以及misclass cost，返回一个total cost，和预测的组
mis_cost_cal <- function(proportion, misclass_cost, method = 'cost', level = NULL){
  # method = cost, pred, both
  # 这里有个robustness的问题，需要proportion是带着names（level）的

  # 标准化proportion
  # 这里还真不能标准化，因为小节点的错误率确实小
  # proportion = proportion / sum(proportion)
  cost_tmp = misclass_cost %*% proportion
  pred = level[which.min(cost_tmp)]
  cost = min(cost_tmp)
  if(method == 'cost'){
    return(cost)
  }else if(method == 'pred'){
    return(pred)
  }else if(method == 'both'){
    return(list(cost = cost, pred = pred))
  }
}

############

nsphere <- function(x){
  # 函数的目的是将K列的矩阵从笛卡尔坐标系变成球坐标系
  # 返回一个K列的矩阵
  # 变换方式：https://www.wikiwand.com/en/N-sphere#/Spherical_coordinates
  if(is.null(dim(x))){
    x <- matrix(x,1)
  }
  K = dim(x)[2]
  if(is.null(K)){
    cat('Something wrong with the input for polar transformation')
  }
  res = x
  res[,1L] = apply(x,1L,function(o_o) sqrt(sum(o_o^2)))
  # Calculate the phi
  phi_func <- function(x_x){
    return(acos(x_x[1] / sqrt(sum(x_x^2))))
  }
  for(i in seq_len(K-1)){
    res[,i+1] = apply(x[,i:K,drop = FALSE],1L,phi_func)
  }
  res[which(x[,K]<0),K] = 2 * pi - res[which(x[,K]<0),K]
  return(res)
}

# 自定义infix函数
`%||%` <- function(lhs, rhs) {
  if (!is.null(lhs)) {
    lhs
  } else {
    rhs
  }
}










