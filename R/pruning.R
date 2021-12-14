#' Title
#'
#' @param fit
#'
#' @return
#' @export
#'
#' @examples
traverse <- function(fit){
  node_saved = fit$treenode
  id_tmp = length(node_saved)
  while(id_tmp >= 1){
    if(!is.null(node_saved[[id_tmp]])){
      node_tmp = node_saved[[id_tmp]][[1]]
      if(!is.na(node_tmp$split_idx)){
        # 中间节点
        node_tmp$leaves = c(node_tmp$idx*2, node_tmp$idx*2 + 1,
                            node_saved[[node_tmp$idx*2]][[1]]$leaves,
                            node_saved[[node_tmp$idx*2 + 1]][[1]]$leaves)
        node_tmp$resub_error = node_saved[[node_tmp$idx*2]][[1]]$resub_error+
          node_saved[[node_tmp$idx*2 + 1]][[1]]$resub_error
        node_tmp$alpha = (node_tmp$misclass - node_tmp$resub_error) /
          (length(node_tmp$leaves) - 1)
        # alpha has to be positive number
        node_tmp$alpha = max(0, node_tmp$alpha,node_saved[[node_tmp$idx*2]][[1]]$alpha,
                             node_saved[[node_tmp$idx*2 + 1]][[1]]$alpha, na.rm = T) # monotonicity
      }else{
        node_tmp$resub_error = node_tmp$misclass
      }
      node_saved[[id_tmp]] = list(node_tmp)
    }
    # print(id_tmp)
    id_tmp = id_tmp - 1
  }
  fit$treenode = node_saved
  return(fit)
}

get_alpha <- function(fit){
  node_saved = fit$treenode
  node_max = length(node_saved)
  alpha_record = rep(NA,node_max)
  for(i in 1:node_max){
    if(!is.null(node_saved[[i]])){
      alpha_record[i] = node_saved[[i]][[1]]$alpha
    }
  }
  return(alpha_record)
}

get_cut_alpha <- function(fit){
  alpha_record <- get_alpha(fit)
  ans_tmp = sort(unique(alpha_record)) # NA被自动干掉了
  if(length(ans_tmp) == 1){ # 说明此时要切掉根结点了
    return(ans_tmp) # 告诉LDATree切的只剩根结点
  }else{
    return(exp(mean(log(ans_tmp[1:2])))) # Geometry average
  }
}

cut_alpha <- function(fit, alpha_tmp){
  node_saved = fit$treenode
  alpha_record <- get_alpha(fit)
  list_tmp = rev(which(alpha_record <= alpha_tmp))
  if(length(list_tmp)!=0){ # 有可能某个子树不需要剪枝
    for(i in 1:length(list_tmp)){
      node_tmp = node_saved[[list_tmp[i]]][[1]]
      # 开始重置为叶子结点
      node_saved[node_tmp$left] <- list(NULL) # set NULL 且不删除
      node_saved[node_tmp$right] <- list(NULL)
      node_tmp$left = node_tmp$right = node_tmp$alpha = node_tmp$criteria =
        node_tmp$split_idx = node_tmp$split_cri = node_tmp$split_na_action = NA
      node_tmp['leaves'] = list(NULL)
      node_saved[[list_tmp[i]]] = list(node_tmp)
    }
    fit$treenode = node_saved
  }
  return(fit)
}

get_mean_se <- function(cv_fit,idx_CV,response,dat,cv_number){
  error_record = numeric(cv_number)
  for(i in 1:cv_number){
    cat('CV:',i,'\n')
    r_tmp = which(idx_CV == i) # 找出第i组的行index
    # 这里要加as.character() 否则apply会把部分factor变成数字，但归根结底还是predict函数不能支持多个x一起，待修改
    predict_tmp = character(length(r_tmp))
    for(j in 1:length(r_tmp)){
      predict_tmp[j] = as.character(predict_LT(cv_fit[[i]],x_new = dat[r_tmp[j],]))
    }
    observe_tmp = response[r_tmp]
    error_record[i] = mean(observe_tmp!=predict_tmp)
  }
  # print(error_record)
  return(c(mean(error_record), sd(error_record) / sqrt(cv_number))) # 开始这里没除sqrt{n} 很愚蠢
}





