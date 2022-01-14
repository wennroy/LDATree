#' Title
#'
#' @param fit
#' @param x_new
#'
#' @return
#' @export
#'
#' @examples
predict_LT <- function(fit, x_new){
  # 别假设了，先把test变量的位置变成和原数据一模一样
  cname_save = fit$cnames
  if(is.null(dim(x_new))){ # 如果传进来的东西是个vector
    x_new = data.frame(matrix(x_new,1))
  }

  # 开始对名字
  matching_position = match(cname_save, colnames(x_new))
  if(anyNA(matching_position)){
    if(dim(x_new)[2] == length(cname_save)){
      colnames(x_new) = cname_save # 假设test的变量位置和原数据一模一样，老罗也是这么做的
    }else{
      print('The new data is missing columns')
      return()
    }
  }else{
    x_new = x_new[,matching_position]
  }

  # 主程序正式启动
  tmp_node = 1
  # 先找到一个最像自己的伙伴, Lazy启动
  friend_list = NULL
  ### 一个距离算法
  while(1){
    # print(tmp_node)
    node_tmp = fit$treenode[[tmp_node]][[1]]
    # print(node_tmp$idx)
    # if(is.null(node_tmp$leaves)){
    if(any(is.na(node_tmp$children))){
      # 到了叶子结点
      if(node_tmp$pred_method == 'mode'){
        # 用众数估计
        return(node_tmp$lda_pred)
      }else{
        # fit LDA
        if(any(is.na(x_new))){
          # 如果有NA，用最像的老哥填充
          idx_NA = which(is.na(x_new))
          if(is.null(friend_list)){ # Lazy 生成 friend_list
            friend_list = best_friend(x_new,fit$dat)
          }
          for(i in idx_NA){
            cursor = 1
            while(1){
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
        return(predict(node_tmp$lda_pred, newdata = x_new)$class)
      }
    }else{
      # 在中间节点，决定下一步往哪里走
      if(pmatch(fit$split_method, c('univariate', 'linear')) == 1){
        inline_x = x_new[1,node_tmp$split_idx] # 即将等待划分的X
      }else{
        cat('Node Index:', node_tmp$idx, '\n')
        # print(node_tmp$linear_split_trans)
        cat('idx_c using:', node_tmp$idx_c, '\n')
        inline_x = node_tmp$linear_split_trans(x_new[1,node_tmp$idx_c])[,node_tmp$split_idx]
      }
      if(class(node_tmp$split_cri) %in% c('numeric', 'integer')){
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
          # 是new level的话也会被分到右面
          tmp_node = 2*tmp_node + ifelse(inline_x %in% node_tmp$split_cri, 0, 1)
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
            tmp_node = 2*tmp_node + ifelse(plugin_tmp %in% node_tmp$split_cri, 0, 1)
          }
        }
      }
    }
  }
}
