#' Splitting in the LDAtree
#'
#' @param x
#' @param y
#' @param mis_curr
#'
#' @return
#' @export
#' @importFrom MASS lda
#' @import dplyr
#'
#' @examples
split_cat <- function(x,y, datx, mis_curr, prior){
  if(!anyNA(x)){
    return(list(split_cat_helper(x,y, datx, mis_curr, prior),NA))
  }else{
    levels(x) = c(levels(x), 'MissingDalian')
    x[is.na(x)] <- 'MissingDalian'
    ans = split_cat_helper(x,y, datx, mis_curr, prior)
    if('MissingDalian' %in% ans){
      return(list(setdiff(ans,'MissingDalian'),1))
    }else{
      return(list(ans,0))
    }
  }
}
split_cat_helper <- function(x,y, datx, mis_curr, prior){
  # Following Loh's 09 paper Procedure A.1

  # x = factor(x) # 这个relevel很有必要，要不然会因为model.matrix有一列全是0从而导致constant
  # 没啥必要，因为如果有一个level是0的话，只会出一个warning
  # Warning message:
  #   In lda.default(x, grouping, ...) : group ADgainai is empty

  ##### 第一种情况 J = 2 #####

  if(length(unique(y)) == 2){ # J = 2, y只有两类
    # cat('Split Cat Situation:1 \n')
    y_class1 = (y == levels(y)[1])
    dat = data.frame(x,y_class1)
    proportion_table = dat %>%
      group_by(x) %>%
      summarise(proportion_r = sum(y_class1))
    proportion_table = proportion_table[order(proportion_table$proportion_r),]

    ans = rep(Inf,nrow(proportion_table)-1) # 记录下错误的个数
    row_tmp = c() # 记录下当前往左节点划分的row有哪些
    for(i in 1:(nrow(proportion_table)-1)){
      row_tmp = c(row_tmp, which(x == proportion_table$x[i]))
      # idx = which(x<= threshold[i])
      y_l = y[row_tmp]
      y_r = y[-row_tmp]
      mis_l = mis_r = 0
      if(length(unique(y_l))!=1){ # 我们采用小于等于的准则
        x_l = datx[row_tmp,]
        # fit = quietly(total_LDA)(x_l, y_l)$result # 用子节点的LDA结果作为划分的依据
        # mis_l = sum(predict(fit,cbind(x_l,y_l))$class != y_l)
        # mis_l = get_error_LDA(x_l, y_l, prior)
        mis_l = pred_LDA(x_l, y_l, prior)[[3]]
      }
      if(length(unique(y_r))!=1){ # 我们采用小于等于的准则
        x_r = datx[-row_tmp,]
        # fit = quietly(total_LDA)(x_r, y_r)$result
        # mis_r = sum(predict(fit,cbind(x_r,y_r))$class != y_r)
        # mis_r = get_error_LDA(x_r, y_r, prior)
        mis_r = pred_LDA(x_r, y_r, prior)[[3]]
      }
      ans[i] = mis_l + mis_r
    }
    idx_threshold = which(ans == min(ans))[1]
    return(proportion_table$x[1:idx_threshold])
  }else if(length(unique(x)) <= 11){
    ##### 第二种情况 n <= 11 #####
    # cat('Split Cat Situation:2 \n')
    level_record = unique(x)
    ans = rep(Inf,2^(length(level_record)-1)-1) # 记录下错误的个数
    row_record = lapply(level_record, function(t) which(x == t)) # 记录下每个level对应的行数
    for(i in 1:length(ans)){
      idx_combn = as.integer(intToBits(i))[1:(length(level_record)-1)]
      row_tmp = unlist(row_record[which(idx_combn==1)])
      # idx = which(x<= threshold[i])
      y_l = y[row_tmp]
      y_r = y[-row_tmp]
      mis_l = mis_r = 0
      if(length(unique(y_l))!=1){ # 我们采用小于等于的准则
        x_l = datx[row_tmp,]
        # fit = quietly(total_LDA)(x_l, y_l)$result # 用子节点的LDA结果作为划分的依据
        # mis_l = sum(predict(fit,cbind(x_l,y_l))$class != y_l)
        # mis_l = get_error_LDA(x_l, y_l, prior)
        mis_l = pred_LDA(x_l, y_l, prior)[[3]]
      }
      if(length(unique(y_r))!=1){ # 我们采用小于等于的准则
        x_r = datx[-row_tmp,]
        # fit = quietly(total_LDA)(x_r, y_r)$result
        # mis_r = sum(predict(fit,cbind(x_r,y_r))$class != y_r)
        # mis_r = get_error_LDA(x_r, y_r, prior)
        mis_r = pred_LDA(x_r, y_r, prior)[[3]]
      }
      ans[i] = mis_l + mis_r
    }
    idx_threshold = which(ans == min(ans))[1]
    return(level_record[which(as.integer(intToBits(idx_threshold))[1:(length(level_record)-1)]==1)])
  }else if((length(unique(x)) > 11) & (length(unique(y)) <= 11)){
    ##### 第三种情况 J <= 11 and n > 20 #####
    # cat('Split Cat Situation:3 \n')
    # Generate X'
    getmode <- function(v){
      uniqv <- unique(v)
      mode_user = uniqv[which.max(tabulate(match(v, uniqv)))]
      return(mode_user)
    }

    dat = data.frame(x,y)
    transfer_matrix = dat %>%
      group_by(x) %>%
      summarise(x_prime = getmode(y))

    x_prime = left_join(data.frame(x), transfer_matrix,by = 'x')[,2]

    level_record = unique(x_prime)
    ans = rep(Inf,2^(length(level_record)-1)-1) # 记录下错误的个数
    row_record = lapply(level_record, function(t) which(x_prime == t)) # 记录下每个level对应的行数
    for(i in 1:length(ans)){
      idx_combn = as.integer(intToBits(i))[1:(length(level_record)-1)]
      row_tmp = unlist(row_record[which(idx_combn==1)])
      # idx = which(x<= threshold[i])
      y_l = y[row_tmp]
      y_r = y[-row_tmp]
      mis_l = mis_r = 0
      if(length(unique(y_l))!=1){ # 我们采用小于等于的准则
        x_l = datx[row_tmp,]
        # fit = quietly(total_LDA)(x_l, y_l)$result # 用子节点的LDA结果作为划分的依据
        # mis_l = sum(predict(fit,cbind(x_l,y_l))$class != y_l)
        # mis_l = get_error_LDA(x_l, y_l, prior)
        mis_l = pred_LDA(x_l, y_l, prior)[[3]]
      }
      if(length(unique(y_r))!=1){ # 我们采用小于等于的准则
        x_r = datx[-row_tmp,]
        # fit = quietly(total_LDA)(x_r, y_r)$result
        # mis_r = sum(predict(fit,cbind(x_r,y_r))$class != y_r)
        # mis_r = get_error_LDA(x_r, y_r, prior)
        mis_r = pred_LDA(x_r, y_r, prior)[[3]]
      }
      ans[i] = mis_l + mis_r
    }
    idx_threshold = which(ans == min(ans))[1]
    result_list = level_record[which(as.integer(intToBits(idx_threshold))[1:(length(level_record)-1)]==1)]
    return(left_join(data.frame(x_prime = result_list),transfer_matrix,by = 'x_prime')[,2])
  }else{
    ##### 第四种情况 Else #####
    # cat('Split Cat Situation:4 \n')
    dummy_matrix = model.matrix(y~x-1) # Get the dummy matrix
    fit = eigen(cov(dummy_matrix)) # Eigen decomposition, 为了防止LDA矩阵不可逆
    eigen_keep = which(round(fit$values,8) > 0) # 保留正值
    X_dummy = dummy_matrix %*% fit$vectors[,eigen_keep] # Projection
    X_dummy = X_dummy[,apply(X_dummy,2,function(x_x) !within_check(y,x_x)), drop = FALSE] # 改掉group constant
    new_data = data.frame(y, X_dummy)
    fit_lda = lda(y~., data = new_data, prior = prior) # 这个LDA需要后面的scaling，所以先不变robust了
    X_num = X_dummy %*% fit_lda$scaling[,1] # Project 到 LD1 上面去
    reexp = unique(data.frame(x,X_num)) # 找到x和X_num的对照表
    threshold = split_noncat(X_num,y,datx,Inf, prior)
    return(reexp$x[which(reexp$X_num <= threshold)])
    }

  # fit = MASS::lda(y~x-1, prior = prior)
  # l = levels(x)[which(coefficients(fit)>0)] # 无法判断是否有提升
  # return(l)
}


# Split_noncat ------------------------------------------------------------

split_noncat <- function(x,y,datx, mis_curr, prior){
  # cat('Split NonCat \n')
  threshold = sort(unique(x))
  # cat('length of x', length(threshold),'\n')
  # ans = ifelse(length(threshold) <= 1000,
  #              split_noncat_small(x,y,datx, mis_curr, prior),
  #              split_noncat_large(x,y,datx, mis_curr, prior))
  if(length(threshold) <= 1000){
    if(!anyNA(x)){
      return(c(split_noncat_small(x,y,datx, mis_curr, prior)[1], NA))
    }else{ # NA 必分到左或者右，所以干脆变成最左和最右
      x1 = x2 = x
      x1[which(is.na(x1))] = min(x, na.rm = T) - 1
      x2[which(is.na(x2))] = max(x, na.rm = T) + 1
      res1 = split_noncat_small(x1,y,datx, mis_curr, prior)
      res2 = split_noncat_small(x2,y,datx, mis_curr, prior)
      return(ifelse(rep(res1[2] <= res2[2],2), c(res1[1],1), c(res2[1],0)))
    }
  }else{
    if(!anyNA(x)){
      return(c(split_noncat_large(x,y,datx, mis_curr, prior)[1], NA))
    }else{ # NA 必分到左或者右，所以干脆变成最左和最右
      x1 = x2 = x
      x1[which(is.na(x1))] = min(x, na.rm = T) - 1
      x2[which(is.na(x2))] = max(x, na.rm = T) + 1
      res1 = split_noncat_large(x1,y,datx, mis_curr, prior)
      res2 = split_noncat_large(x2,y,datx, mis_curr, prior)
      return(ifelse(rep(res1[2] <= res2[2],2), c(res1[1],1), c(res2[1],0)))
    }
  }
  # return(ans)
}

split_noncat_small <- function(x,y,datx, mis_curr, prior){ # 这一步跑得太太太太慢
  # cat('Split NonCat Small \n')
  threshold = sort(unique(x)) # sort会自动干掉NA
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
      # mis_l = get_error_LDA(x_l, y_l, prior)
      mis_l = pred_LDA(x_l, y_l, prior)[[3]]
    }
    if(length(unique(y_r))!=1){ # 我们采用小于等于的准则
      x_r = datx[-idx,]
      # fit = quietly(total_LDA)(x_r, y_r)$result
      # mis_r = sum(predict(fit,cbind(x_r,y_r))$class != y_r)
      # mis_r = get_error_LDA(x_r, y_r, prior)
      mis_r = pred_LDA(x_r, y_r, prior)[[3]]
    }
    ans[i] = mis_l + mis_r
  }
  idx_threshold = which(ans == min(ans))[1]
  return(c(threshold[idx_threshold], min(ans)))
  # cat(ans,mis_curr,'\n')
  # print(c(ans[idx_threshold],mis_curr))
  # if(ans[idx_threshold] >= mis_curr){
  #   return(NULL)
  # }else{
  #   return(threshold[idx_threshold])
  # }
}

split_noncat_large <- function(x,y,datx, mis_curr, prior){ # 这一步跑得不够稳健
  # cat('Split NonCat Large\n')
  threshold = sort(unique(x))
  left_pointer = 1
  right_pointer = length(threshold)
  ans = matrix(c(0,Inf),1,2)
  while(right_pointer >= left_pointer){
    # cat(left_pointer, right_pointer, '\n')
    current_index = (left_pointer + right_pointer) %/% 2
    idx = which(x<= threshold[current_index])
    y_l = y[idx]
    y_r = y[-idx]
    mis_l = mis_r = 0
    if(length(unique(y_l))!=1){ # 我们采用小于等于的准则
      x_l = datx[idx,]
      # fit = quietly(total_LDA)(x_l, y_l)$result # 用子节点的LDA结果作为划分的依据
      # mis_l = sum(predict(fit,cbind(x_l,y_l))$class != y_l)
      # mis_l = get_error_LDA(x_l, y_l, prior)
      mis_l = pred_LDA(x_l, y_l, prior)[[3]]
    }
    if(length(unique(y_r))!=1){ # 我们采用小于等于的准则
      x_r = datx[-idx,]
      # fit = quietly(total_LDA)(x_r, y_r)$result
      # mis_r = sum(predict(fit,cbind(x_r,y_r))$class != y_r)
      # mis_r = get_error_LDA(x_r, y_r, prior)
      mis_r = pred_LDA(x_r, y_r, prior)[[3]]
    }
    ans = rbind(ans, c(current_index, mis_l + mis_r))
    if(mis_l >= mis_r){
      right_pointer = current_index - 1
    }else{
      left_pointer = current_index + 1
    }
  }
  idx_threshold = ans[which.min(ans[,2]),1]
  # print(ans)
  return(c(threshold[idx_threshold], min(ans[,2])))
  # if(ans[idx_threshold] >= mis_curr){
  #   return(NULL)
  # }else{
  #   return(threshold[idx_threshold])
  # }
}


split_fact_uni <- function(x,y, prior){
  # With given prior
  # fit_split = lda(y~x)
  # print(fit_split)
  # cat(fit_split$means,prior,'\n')

  # gm = data.frame(y,x) %>%
  #   group_by(y) %>%
  #   summarise(x = mean(x))

  gm = sapply(levels(y),function(o_o) mean(x[y==o_o], na.rm = TRUE))
  # print(fit_split$means)
  prior = prior[prior!=0] # 去掉那些空的组
  # gm_obs = cbind(fit_split$means,prior)
  gm_obs = cbind(gm[!is.nan(gm)],prior) # 那些空组就去掉了
  # gm_obs = cbind(gm[,2],prior)
  gm_obs = gm_obs[order(gm_obs[,1]),] # 从小到大排序
  # print(gm_obs)
  # weighted average: (nj-1)sigma2
  sigma2_tmp = weighted.mean(tapply(x,y,FUN = var),w = sapply(table(y),FUN = function(o_o) max(o_o-1,0)))


  # 先找到choose k 2个交点
  possible_cut = c()
  for(o_o in 1:(dim(gm_obs)[1]-1)){
    for(o_0 in (o_o+1):dim(gm_obs)[1]){
      possible_cut = c(possible_cut, (gm_obs[o_o,1] + gm_obs[o_0,1]) / 2 + sigma2_tmp * log(gm_obs[o_o,2] / gm_obs[o_0,2]) / (gm_obs[o_0,1] - gm_obs[o_o,1]))
    }
  }
  # 下面是每个LDA的cut
  # print(possible_cut)
  candidate_cut = sort(unique(possible_cut)[range(x)[1] < unique(possible_cut) & unique(possible_cut) < range(x)[2]])
  if(length(candidate_cut)==0){ # 全部的分割点都在外面，停止分割
    return(NULL)
  }
  plug_in_value = sort(c(range(x),candidate_cut))[-length(candidate_cut)-2] + diff(sort(c(range(x),candidate_cut)))/2
  # print(plug_in_value)
  res_tmp = sapply(1:dim(gm_obs)[1], function(o_o) (plug_in_value - gm_obs[o_o,1])^2 - 2 * sigma2_tmp * log(gm_obs[o_o,2]))
  # print(res_tmp)
  res_tmp = apply(res_tmp,1,which.min)

  final_cut = c() # 真正的cut
  for(o_o in 2 : length(res_tmp)){
    if(res_tmp[o_o] != res_tmp[o_o-1]){
      final_cut = c(final_cut, candidate_cut[o_o - 1])
    }
  }
  if(length(final_cut) == 0){
    return(NULL)
  }
  # 最后检查一下有没有空组
  final_cut_pro = c()
  test_tmp = table(cut(x,breaks = c(-Inf,final_cut,Inf)))
  # print(test_tmp)
  for(o_o in 1 : length(final_cut)){
    if(test_tmp[o_o] != 0){
      final_cut_pro = c(final_cut_pro, final_cut[o_o])
    }
  }
  if(test_tmp[length(test_tmp)] == 0){
    final_cut_pro = final_cut_pro[-length(final_cut_pro)]
  }
  if(length(final_cut_pro) == 0){
    return(NULL)
  }
  return(final_cut_pro)
}







