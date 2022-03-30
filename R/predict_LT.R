#' Title
#'
#' @param fit
#' @param x_new
#'
#' @return
#' @export
#'
#' @examples
predict_LT <- function(fit, x_new, dat = NULL){
  # 这个function只支持一行的输入
  # 如果要实现多行的输入，要使用apply系列的函数

  # 传递参数dat主要是为了NA的情况
  if(any(is.na(x_new)) & is.null(dat)){
    stop('For prediction with missing value, please provide the original data for imputation')
  }

  # 别假设了，先把test变量的位置变成和原数据一模一样
  cname_save = fit$cnames
  if(is.null(dim(x_new))){ # 如果传进来的东西是个vector
    x_new = data.frame(matrix(x_new,1))
  }

  # 先将传入的dat变成只有X的版本
  if(!is.null(dat)){
    if(is.null(colnames(dat)) & (dim(dat)[2] != length(cname_save))){
      stop('The data entry has something wrong. Column name is missing, and we do not know which coliumn represents the response')
    }
    matching_position = match(cname_save, colnames(dat))
    if(anyNA(matching_position)){
      stop('The data entry has something wrong. Some columns might be missing.')
    }
    dat = dat[,matching_position]
  }

  # 开始对名字
  matching_position = match(cname_save, colnames(x_new))
  if(anyNA(matching_position)){
    if(dim(x_new)[2] == length(cname_save)){
      colnames(x_new) = cname_save # 假设test的变量位置和原数据一模一样，老罗也是这么做的
    }else{
      stop('The new data has the wrong number of columns')
    }
  }else{
    x_new = x_new[,matching_position]
  }

  # 主程序正式启动
  tmp_node = 1
  # 先找到一个最像自己的伙伴, Lazy启动
  friend_list = NULL
  ### 一个距离算法

  cov_class_new = sapply(x_new,class) %in% c('numeric', 'integer') # 为了后面检测new level and NA
  ##### 检查new_level和NA #####
  # FACT的话要将所有的x_new都变成数字
  if(fit$select_method == 'FACT'){
    # cov_class_new = sapply(x_new,class) %in% c('numeric', 'integer')
    # cat('cov_class:',cov_class_new,'\n')
    if(!all(cov_class_new)){
      # 如果不全是数字
      for(o_o in which(!cov_class_new)){
        cat_trans_tmp = fit$cat_trans[[fit$cnames[o_o]]]
        if(unlist(x_new[o_o]) %in% cat_trans_tmp[,1]){
          # 既不是NA也不是new_level
          # 下面这一行有as.character 因为如果两个factor level不同的时候
          # 会报错误 level sets of factors are different
          x_new[o_o] = cat_trans_tmp[cat_trans_tmp[,1] == as.character(unlist(x_new[o_o])),2]
        }else{
          x_new[o_o] = NA # 这可能会导致出现NA，给后面带来麻烦
        }
      }
    }

    # 现在从FACT出来的应该都是数字+NA了
    # 这里产生一个分水岭，对于能Handle NA的方法，不进行任何操作
    # 对于无法Handle NA的方法，这里直接填上去。
    if(any(is.na(x_new))){
      # print(x_new)
      x_new = class_centroid_impute(dat, fit$response, prior = fit$prior, cov_class = fit$cov_class, cat_trans = fit$cat_trans,
                                    type = 'all', x_new = x_new)
    }
  }

  # LDA的话，New level -> NA
  # NA -> 先随机变成Mean/Mode好了
  # 先检测new level
  if(fit$select_method == 'LDATree'){
    if(!all(cov_class_new)){
      for(o_o in which(!cov_class_new)){
        # print(x_new)
        # print(o_o)
        if(!unlist(x_new[o_o]) %in% names(fit$level_record[[o_o]])){
          # 不加unlist的话会报莫名其妙的错误
          # 让一个1*1的df Female %in% c('Female') = FALSE
          x_new[o_o] = names(fit$level_record[[o_o]])[which.max(fit$level_record[[o_o]])]
        }
      }
    }
    for(o_o in seq(x_new)){
      if(is.na(x_new[o_o])){
        # Numerical
        if(class(fit$level_record[[o_o]]) == 'numeric'){
          x_new[o_o] = fit$level_record[[o_o]]
        }else{
          x_new[o_o] = names(fit$level_record[[o_o]])[which.max(fit$level_record[[o_o]])]
        }
      }
    }
  }
  # # 再检测是否有NA
  # 但这会破坏原有的结构，先放在这里好了
  # 因为之前对于NA的数据我们也是可以handle的
  # 2022/02/23 但是现在我还是准备先忽略这个问题
  # 因为日后还是要修改这个NA imputing

  repeat{
    node_tmp = fit$treenode[[tmp_node]][[1]]
    # if(is.null(node_tmp$leaves)){
    if(any(is.na(node_tmp$children))){
      # 到了叶子结点
      if(node_tmp$pred_method == 'mode'){
        # 用众数估计
        return(node_tmp$node_pred)
      }else{
        # fit LDA
        # 这个循环内的代码可以变一变
        if(any(is.na(x_new))){
          # 如果有NA，用最像的老哥填充
          idx_NA = which(is.na(x_new))
          if(is.null(friend_list)){ # Lazy 生成 friend_list
            friend_list = best_friend(x_new,fit$dat)
          }
          for(i in idx_NA){
            cursor = 1
            repeat{
              if(!is.na(fit$dat[friend_list[cursor],i])){
                x_new[i] = fit$dat[friend_list[cursor],i]
                break
              }else{
                cursor = cursor + 1
              }
            }
          }
        }
        # 扔进LDA
        return(predict(node_tmp$node_pred, newdata = x_new)$class)
      }
    }else{ # 在中间节点，决定下一步往哪里走
      inline_x = x_new[1,node_tmp$split_idx] # 即将等待划分的X
      dim(inline_x) <- NULL # 因经常弄出来一个矩阵，导致后面的代码跑不动
      if(pmatch(fit$split_method, c('univariate', 'linear')) == 2){
        inline_x = node_tmp$linear_split_trans(x_new[1,node_tmp$idx_c])[,node_tmp$split_idx]
      }
      if(fit$cov_class[node_tmp$split_idx]){
        # continous variable
        if(!is.na(inline_x)){
          # tmp_node = 2*tmp_node + ifelse(inline_x <= node_tmp$split_cri, 0, 1)
          # 下面开始写的是大于等于，但是发现不太对劲，应该是大于，因为要找小于等于的反面
          tmp_node = node_tmp$children[which.min(inline_x > c(node_tmp$split_cri,Inf))]
        }else{
          if(!is.na(node_tmp$split_na_action)){
            # 如果划分的时候包含了NA的情况
            tmp_node = 2*tmp_node + ifelse(node_tmp$split_na_action, 0, 1)
          }else{
            # 如果没用NA划分，但是出现了NA
            if(is.null(friend_list)){ # Lazy 生成 friend_list
              friend_list = best_friend(x_new,fit$dat)
            }
            # 找到第一个不是NA的plugin
            plugin_tmp = NULL
            cursor = 1
            while(is.null(plugin_tmp)){
              if(!is.na(fit$dat[friend_list[cursor],node_tmp$split_idx])){
                plugin_tmp = fit$dat[friend_list[cursor],node_tmp$split_idx]
                break
              }else{
                cursor = cursor + 1
              }
            }
            tmp_node = 2*tmp_node + ifelse(plugin_tmp <= node_tmp$split_cri, 0, 1)
          }
        }
      }else{
        # Categorical Variable
        if(!is.na(inline_x)){
          # new level的话先被打成NA，然后再操作
          # tmp_node = 2*tmp_node + ifelse(inline_x %in% node_tmp$split_cri, 0, 1)
          if(fit$select_method == 'FACT'){
            # 对于FACT来说，我们在过程中一直保持numerical，直到最后画图再变回去。
            # FACT 对于自己内部CV剪枝的时候，cat也还是numeric，只有在预测新数据的时候，
            # FACT才需要做cat到num的转换
            if(!(class(inline_x) %in% c('numeric', 'integer'))){
              cat_trans_tmp = fit$cat_trans[[fit$cnames[node_tmp$split_idx]]]
              inline_x = cat_trans_tmp[cat_trans_tmp[,1] == inline_x,2]
            }
            tmp_node = node_tmp$children[which.min(inline_x > c(node_tmp$split_cri,Inf))]
            # tmp_node = 2*tmp_node + ifelse(inline_x %in% node_tmp$split_cri, 0, 1)
          }else if(fit$select_method == 'LDATree'){
            # For LDATree, the criteria is a vector of levels
            tmp_node = 2*tmp_node + ifelse(inline_x %in% node_tmp$split_cri, 0, 1)
          }else{
            stop('Prediction for other methods is still under development')
          }
        }else{ # 如果是NA的话
          if(!is.na(node_tmp$split_na_action)){
            # 如果划分的时候包含了NA的情况
            tmp_node = 2*tmp_node + ifelse(node_tmp$split_na_action, 0, 1)
          }else{
            # 如果没用NA划分，但是出现了NA
            if(is.null(friend_list)){ # Lazy 生成 friend_list
              friend_list = best_friend(x_new,fit$dat)
            }
            # 找到第一个不是NA的plugin
            plugin_tmp = NULL
            cursor = 1
            while(is.null(plugin_tmp)){
              if(!is.na(fit$dat[friend_list[cursor],node_tmp$split_idx])){
                plugin_tmp = fit$dat[friend_list[cursor],node_tmp$split_idx]
                break
              }else{
                cursor = cursor + 1
              }
            }
            tmp_node = 2*tmp_node + ifelse(plugin_tmp %in% node_tmp$split_cri, 0, 1)
          }
        }
      }
    }
  }
}

predict_Tree <- function(fit, x_new, data = NULL){
  # wrapper for prediction
  if(is.null(dim(x_new)) | dim(x_new)[1] == 1){
    return(predict_LT(fit, x_new, data))
  }else{
    # Classification Tree, 结果一定是factor
    res = character(dim(x_new)[1])
    for(i in seq(res)){
      res[i] = as.character(predict_LT(fit, x_new[i,], data))
    }
    return(res)
  }
}






