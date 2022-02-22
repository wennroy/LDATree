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
# These functions are used for NA imputation

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

class_centroid_impute <- function(xs, y, prior, cov_class = NULL, cat_trans = NULL,
                                  type = 'all', x_new = NULL){
  # 如果传进来的xs有categorical的变量
  # 要先变成numerical: 这一步的目的是只改变xs
  # 这个在predict的时候有时候会遇到
  # 尤其是CV的时候，要同时predict很多行
  # 2022/02/11 停在这里
  xs_check = sapply(xs,class)
  if('factor' %in% xs_check){
    idx_trans = which(xs_check == 'factor')
    for(i in idx_trans){
      trans_table = cat_trans[[colnames(xs)[i]]]
      xs[,i] = left_join(data.frame(x = xs[,i]),trans_table, by = "x")[,2]
    }
  }

  # type = c('all', 'single')
  # 这个function以后可以变成一个很通用的function，专治各种missing问题
  # 这个函数有点慢
  print('Imputing data...')
  cov_class <- cov_class %||% (sapply(xs,class) %in% c('numeric', 'integer')) # 如果没给，那就用数据本来的样子
  if(is.null(dim(xs))){# 传进来的如果是一个vector，把它变成matrix，才可以使用apply
    xs = matrix(xs, nrow = length(xs), ncol = 1)
  }

  # 每一列是一个类别
  m_table <- sapply(levels(y), function(x_x) apply(xs,2,function(o_o) mean(o_o[y == x_x], na.rm = TRUE)))
  v_table <- sapply(levels(y), function(x_x) apply(xs,2,function(o_o) var(o_o[y == x_x], na.rm = TRUE)))
  print('I am running')
  # 如果方差为0，则把方差变成一个比较小的数字
  v_table[v_table < 1e-5] = 1e-5
  # 对于categorical变量，因为数值是离散的，
  # 所以把mean变成距离最近的类别的数值
  if(!all(cov_class)){
    for(i in which(!cov_class)){
      cat_trans_tmp = cat_trans[[names(xs)[i]]]
      for(j in 1:ncol(m_table)){
        # print(cat_trans_tmp[which.min(abs(m_table[i,j] - cat_trans_tmp[,2])),2])
        m_table[i,j] = cat_trans_tmp[which.min(abs(m_table[i,j] - cat_trans_tmp[,2])),2]
      }
    }
  }

  # 传入一个x，计算K个中心，然后给出估计
  normal_density_impute <- function(x, m_table, v_table, prior){
    if(all(is.na(x))){ # 如果这个数据所有列都missing，随便选一个centroid给他塞进去好了。
      return(m_table[,sample(ncol(m_table),1)])
    }
    # 选出那些没有Missing的位置
    idx_missing = which(is.na(x))
    # x_now = x[-idx_missing]
    # m_table_now = m_table[-idx_missing,]
    # v_table_now = v_table[-idx_missing,]
    dist_saved = sapply(seq(ncol(m_table)), function(o_o) -log(prior[o_o]) + 0.5 * sum(log(v_table[-idx_missing,o_o]) +(x[-idx_missing] - m_table[-idx_missing,o_o])^2 / v_table[-idx_missing,o_o]))
    # if(any(is.nan(dist_saved))){
    #   print(m_table)
    #   print(v_table)
    #   print(prior)
    #   print(idx_missing)
    #   print(dist_saved)
    #   stop('I am here!')
    # }
    x[idx_missing] = m_table[idx_missing,which.min(dist_saved)]
    return(x)
  }


  # 只预测一个新的X
  if(type == 'single'){
    x_new = normal_density_impute(x_new, m_table, v_table, prior)
    return(x_new)
  }

  # 对数据的每一行进行填充
  for(i in 1:nrow(xs)){
    if(any(is.na(xs[i,]))){
      xs[i,] = normal_density_impute(xs[i,], m_table, v_table, prior)
    }
  }

  return(xs)
}


