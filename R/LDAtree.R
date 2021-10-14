#' Main function in LDATree
#'
#' Fit the LDATree, and return the treenodes.
#'
#' @param formula an object of class 'formula'
#' @param data Data frame with response variable in it.
#'
#' @return List
#' @export
#'
#' @examples
LDAtree <- function(formula, data, prior = NULL, max_level = 5, min_nsize = NULL){
  # prior 的顺序需要和data中的levels相同
  # data 是含有response的， dat不含有
  response = model.frame(formula, data)[,1L]
  if(!is.factor(response)){
    response = as.factor(response)
  }

  # 去掉y是NA的那些数据
  idx_notNA = which(!is.na(response))
  response = response[idx_notNA]
  data = data[idx_notNA,]

  if(is.null(min_nsize)){
    min_nsize = max(length(response) %/% 100, 5) # 每个node最少5个
  }

  # Design Matrix的方法很不对！因为这样卡方选变量根本跑不动了。
  # 我们需要找个办法找到所有涉及的变量。
  # dat = as.data.frame(model.matrix(formula, data)[,-1])
  dat = data[,sapply(labels(terms(formula, data = data)),function(x) which(x == colnames(data)))]
  col_idx = 1:ncol(dat)

  # Put prior in somewhere
  if(is.null(prior)){
    prior = tabulate(response)/nrow(data) # 如果没给prior，用estimated prior
  }

  # 下面就是陈年老code

  queue = list(list(1:nrow(dat),col_idx,1L)) # Used to save the current intermediate nodes
  c_name = colnames(dat) # Used to denote the criteria
  node_saved = list() # Save the nodes for printing the tree

  while(length(queue)!=0){ # 当还有节点没有被划分好
    node_tmp = generate_node(idx_r = queue[[1]][[1]],
                             idx_c = queue[[1]][[2]],
                             idx = queue[[1]][[3]]) # Get the current subset index
    cat('The current node index is', node_tmp$idx, '\n')
    queue = queue[-1] # Remove the idx from the waiting list

    if(node_tmp$size == 0){ # If there is no data in this current node
      # 我现在心生疑惑，感觉这个情况永远不会遇到。
      # numeric的划分应该不会遇到，只可能是LDA分的时候系数都是一个符号？
      # 如果真的有机会跑到这里我再优化吧
      node_tmp$misclass = 0
      node_saved[[node_tmp$idx]] = list(node_tmp)
      next # 出口1
    }

    response_tmp = response[node_tmp$idx_r] # get the subset of the response

    ## count开始，数清楚每一类有多少个
    node_tmp$portion = numeric(length(levels(response)))
    for(i in 1:length(node_tmp$portion)){
      node_tmp$portion[i] = sum(response_tmp == levels(response)[i])
    }
    ## count结束了


    Jt = sum(node_tmp$portion != 0) # number of class in node t 当前节点有几类
    if(Jt == 1 || node_tmp$covs == 0){ # If all of the data belongs to one class, or there is no covariates left
      node_tmp$misclass = node_tmp$size - max(node_tmp$portion) # 既然分不下去了，那就只能预测众数
      node_tmp$pred_method = 'mode'
      node_tmp$lda_pred = levels(response)[which.max(node_tmp$portion)]
      node_saved[[node_tmp$idx]] = list(node_tmp)
      next # 出口2
    }


    # 产生所需要的数据
    dat_tmp = data.frame(dat[node_tmp$idx_r,node_tmp$idx_c]) # Generate temporary data
    # 开始治理那些所有除了NA之外值都一样的列
    idx_c_keep = apply(dat_tmp,2,function(x) length(unique(x[!is.na(x)])) > 1)
    node_tmp$idx_c = node_tmp$idx_c[idx_c_keep]
    node_tmp$covs = length(node_tmp$idx_c)
    dat_tmp = data.frame(dat[node_tmp$idx_r,node_tmp$idx_c]) # Generate temporary data
    colnames(dat_tmp) = c_name[node_tmp$idx_c] # rename the data columns 这一步可能是为了增加robustness？

    if(node_tmp$covs == 0){ # If there is no covariates left
      node_tmp$misclass = node_tmp$size - max(node_tmp$portion) # 既然分不下去了，那就只能预测众数
      node_tmp$pred_method = 'mode'
      node_tmp$lda_pred = levels(response)[which.max(node_tmp$portion)]
      node_saved[[node_tmp$idx]] = list(node_tmp)
      next # 出口3
    }

    # fit <- total_LDA(dat_tmp, response_tmp)
    # predict_tmp = predict(fit,newdata = dat_tmp) # Predict using the fitted model
    # node_tmp$misclass = sum(predict_tmp$class != response_tmp)
    node_tmp$misclass = get_error_LDA(dat_tmp, response_tmp, prior) # 这个信息会一直保留，因为最后需要展示在图上

    # 判断是否到达了 maximum level，或者minimum node size，到达了便退出
    if(node_tmp$idx >= 2 ** (max_level-1)){
      ans = pred_LDA(dat_tmp, response_tmp, prior)
      node_tmp$pred_method = ans[[1]]
      node_tmp$lda_pred = ans[[2]]
      node_saved[[node_tmp$idx]] = list(node_tmp)
      next # 出口4
    }

    if(node_tmp$misclass > 0){ # If there is no perfect seperation, we choose one of the best seperator
      # 能进这里来的，肯定是可以跑LDA function的
      # We should also record the covariates
      # 不管是什么类型的数据，我们都可以复用

      # Variable selection
      chi_stat = apply(dat_tmp,2,function(x) var_select_all(x,response_tmp, node_tmp$size, Jt))
      c_split = which(chi_stat == max(chi_stat))[1] # 这个是新的index，不是旧的


      # Splitting

      flag_class = class(dat_tmp[,c_split]) %in% c('numeric', 'integer')
      if(flag_class){
        threshold = split_noncat(dat_tmp[,c_split],response_tmp,dat_tmp, node_tmp$misclass, prior)
        node_tmp$split_cri = threshold[1]
        node_tmp$split_na_action = threshold[2]
        # if(is.null(threshold)){ # 如果没有提升
        #   # 这里的代码，对于预测新数据很关键
        #   # 暂时我们可以先允许早期停止，之后的话，加上了pruning之后，我们可以选择不停止
        #   ans = pred_LDA(dat_tmp, response_tmp)
        #   node_tmp$pred_method = ans[[1]]
        #   node_tmp$lda_pred = ans[[2]]
        #   node_saved[[node_tmp$idx]] = list(node_tmp)
        #   next # 出口5
        # }else{
        #   idx_left = which(dat_tmp[,c_split] <= threshold)
        #   idx_right = setdiff(1:node_tmp$size,idx_left)
        # }
        idx_left = which(dat_tmp[,c_split] <= threshold[1])
        # 如果有NA且被分到左面的话：
        if(anyNA(dat_tmp[,c_split]) & node_tmp$split_na_action == 1){
          idx_left = c(idx_left, which(is.na(dat_tmp[,c_split])))
        }
        idx_right = setdiff(1:node_tmp$size,idx_left)
      }else{
        left_group = split_cat(dat_tmp[,c_split],response_tmp, dat_tmp, node_tmp$misclass, prior)
        node_tmp$split_cri = left_group[[1]]
        node_tmp$split_na_action = left_group[[2]]
        idx_left = which(dat_tmp[,c_split] %in% left_group[[1]])
        if(anyNA(dat_tmp[,c_split]) & node_tmp$split_na_action == 1){
          idx_left = c(idx_left, which(is.na(dat_tmp[,c_split])))
        }
        idx_right = setdiff(1:node_tmp$size,idx_left)
      }
      subnode_index_c = node_tmp$idx_c
      node_tmp$split_idx = node_tmp$idx_c[c_split]

      # If both sides have data (ortherwise, this node should be a leaf node)
      # 下面这个看起来也是一个永远不会遇到的情况

      if(length(idx_left) != 0 & length(idx_right) != 0){
        # Define the spliting rule for numeric variable
        node_tmp$criteria = ifelse(flag_class, paste(c(c_name[node_tmp$split_idx],'\u2264',
                                                       ifelse(is.na(node_tmp$split_na_action),'',
                                                              ifelse(node_tmp$split_na_action,'**','++')),
                                                       threshold[1]),collapse = ' ')
                                   ,paste(c(c_name[node_tmp$split_idx],'in {', paste(left_group,sep = ', '),
                                            ifelse(is.na(node_tmp$split_na_action),'',
                                                   ifelse(node_tmp$split_na_action,'**','++')), '}'),collapse = ' '))

        node_tmp$left = 2 * node_tmp$idx # Decide the child index for the current node
        queue = append(queue,list(list(node_tmp$idx_r[idx_left],subnode_index_c, node_tmp$left)))
        node_tmp$right = 2 * node_tmp$idx + 1 # Decide the child index for the current node
        queue = append(queue,list(list(node_tmp$idx_r[idx_right],subnode_index_c, node_tmp$right)))
      }
    }
    node_saved[[node_tmp$idx]] = list(node_tmp) # 出口6
  }
  res = list()
  res$formula = formula
  res$treenode = node_saved
  res$dat = dat
  res$prior = prior
  res$response = response # 为了后面的画图
  res$response_name = colnames(model.frame(formula, data))[1]
  cat('The LDA tree is completed.\n')
  return(res)
}
