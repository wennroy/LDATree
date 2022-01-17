fact <- function(formula, data, prior = NULL, max_level = 10,
                    min_nsize = NULL, cv_number = 10){
  # prior 的顺序需要和data中的levels相同
  # data 是含有response的， dat不含有


  # 数据预处理 -------------------------------------------------------------------


  ### 将response变成factor
  response = model.frame(formula, data)[,1L]
  if(!is.factor(response)){
    response = as.factor(response)
  }
  ###

  ### 去掉y是NA的那些数据
  idx_notNA = which(!is.na(response))
  response = response[idx_notNA]
  data = data[idx_notNA,]
  ###

  ### 设定每个节点的最小数据量
  if(is.null(min_nsize)){
    # min_nsize = max(length(response) %/% 100, 5) # 每个node最少5个
    min_nsize = 5 # FACT
  }
  ###

  ### 设定Prior
  if(is.null(prior)){
    prior = tabulate(response)/nrow(data) # 如果没给prior，用estimated prior
  }
  ###

  ### 找到所有要使用的X，合在一起变成dat
  dat = data[,sapply(labels(terms(formula, data = data)),function(x) which(x == colnames(data)))]
  ###


  # Model 部分 ----------------------------------------------------------------
  fit = tree_growing_fact(response, dat, prior, max_level, min_nsize, select_method)

  # Pruning -----------------------------------------------------------------
#
#   T_saved = fit$treenode
#   fit = traverse(fit) # Get the alpha
#   # Divide the data into ten parts
#   set.seed(dim(dat)[1]*dim(dat)[2]) # Fixed seed
#   idx_CV = sample(c(rep(1:(cv_number-1),each = dim(dat)[1]%/%cv_number),
#                     numeric(dim(dat)[1] - (cv_number-1)* (dim(dat)[1]%/%cv_number))+cv_number))
#   cv_fit = vector("list", cv_number)
#   for(i in 1:cv_number){
#     r_tmp = which(idx_CV == i) # 找出第i组的行index
#     cv_fit[[i]] = tree_growing(response[-r_tmp], dat[-r_tmp,], prior, max_level, min_nsize) # 找出母树
#     cv_fit[[i]] = traverse(cv_fit[[i]]) # update alpha
#   }
#   # plotall(cv_fit[[2]])
#   CV_table = data.frame(Tree = numeric(0),
#                         Tnodes = numeric(0),
#                         Mean_MSE = numeric(0),
#                         SE_Mean = numeric(0),
#                         alpha = numeric(0))
#   CV_table[dim(CV_table)[1]+1,] = c(dim(CV_table)[1]+1,sum(!sapply(fit$treenode,is.null)),
#                                     get_mean_se(cv_fit,idx_CV,response,dat,cv_number),-1)
#   while(sum(!sapply(fit$treenode,is.null)) > 1){ # 除非是切到根结点了
#     cat('The current number of node is', sum(!sapply(fit$treenode,is.null)), '\n')
#     # 找到切割的alpha
#     alpha_tmp = get_cut_alpha(fit)
#     fit = cut_alpha(fit,alpha_tmp)
#     fit = traverse(fit)
#     for(i in 1:cv_number){
#       cv_fit[[i]] = cut_alpha(cv_fit[[i]],alpha_tmp)
#       cv_fit[[i]] = traverse(cv_fit[[i]])
#     }
#     # plotall(cv_fit[[3]])
#     CV_table[dim(CV_table)[1]+1,] = c(dim(CV_table)[1]+1,sum(!sapply(fit$treenode,is.null)),
#                                       get_mean_se(cv_fit,idx_CV,response,dat,cv_number),alpha_tmp)
#   }
#   # Evaluation --------------------------------------------------------------
#   CV_table = CV_table[dim(CV_table)[1]:1,] # 反向排序，因为which.min会选择最上面的，
#   # 但我们需要最小的树
#   k_se = 1
#   idx_star = which.min(CV_table$Mean_MSE)
#   star_threshold = CV_table$Mean_MSE[idx_star] + k_se * CV_table$SE_Mean[idx_star]
#   alpha_final = CV_table$alpha[which(CV_table$Mean_MSE<=star_threshold)[1]]
#   fit$treenode = T_saved
#   fit = traverse(fit) # Get the alpha
#   CV_table = CV_table[dim(CV_table)[1]:1,]
#   tmp = 1
#   while(CV_table$alpha[tmp] <= alpha_final){
#     fit = cut_alpha(fit,CV_table$alpha[tmp])
#     fit = traverse(fit)
#     tmp = tmp + 1
#   }
#
#   fit$response_name = colnames(model.frame(formula, data))[1]
#   fit$CV_table = CV_table[dim(CV_table)[1]:1,]
  return(fit)
}



# Tree Growing Function ---------------------------------------------------



tree_growing_fact <- function(response, dat, prior, max_level, min_nsize, select_method){
  # Model 部分 ----------------------------------------------------------------
  col_idx = 1:ncol(dat)
  queue = list(list(1:nrow(dat),col_idx,1L)) # Used to save the current intermediate nodes
  c_name = colnames(dat) # Used to denote the criteria
  node_saved = list() # Save the nodes for printing the tree
  Nj = table(response)
  no_j = length(unique(response))

  while(length(queue)!=0){ # 当还有节点没有被划分好

    ### 拿到当前节点的信息
    node_tmp = generate_node_fact(idx_r = queue[[1]][[1]],
                             idx_c = queue[[1]][[2]],
                             idx = queue[[1]][[3]]) # Get the current subset index
    cat('The current node index is', node_tmp$idx, '\n')
    queue = queue[-1] # Remove the idx from the waiting list
    response_tmp = response[node_tmp$idx_r] # get the subset of the response
    ###

    ### 先拿到一个response的大致分布
    # node_tmp$portion = numeric(length(levels(response)))
    # for(i in 1:length(node_tmp$portion)){
    #   node_tmp$portion[i] = sum(response_tmp == levels(response)[i])
    # }
    node_tmp$portion = as.numeric(table(response_tmp))
    ###

    ### 出口1: 无X，或者所有Y都相同
    Jt = sum(node_tmp$portion != 0) # number of class in node t 当前节点有几类
    if(Jt == 1 | node_tmp$covs == 0){ # If all of the data belongs to one class, or there is no covariates left
      node_tmp$misclass = node_tmp$size - max(node_tmp$portion) # 既然分不下去了，那就只能预测众数
      node_tmp$pred_method = 'mode'
      node_tmp$node_pred = levels(response)[which.max(node_tmp$portion)]
      # 这里要改成prior版本的，带misclassification cost
      node_saved[[node_tmp$idx]] = list(node_tmp)
      next # 出口1
    }
    ###


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
      node_tmp$node_pred = levels(response)[which.max(node_tmp$portion)]
      node_saved[[node_tmp$idx]] = list(node_tmp)
      next # 出口2
    }

    # fit <- total_LDA(dat_tmp, response_tmp)
    # predict_tmp = predict(fit,newdata = dat_tmp) # Predict using the fitted model
    # node_tmp$misclass = sum(predict_tmp$class != response_tmp)

    # node_tmp$misclass = get_error_LDA(dat_tmp, response_tmp, prior) # 这个信息会一直保留，因为最后需要展示在图上
    # ans = pred_LDA(dat_tmp, response_tmp, prior)
    node_tmp$misclass = node_tmp$size - max(node_tmp$portion) # 这个信息会一直保留，因为最后需要展示在图上
    # node_tmp$pred_method = ans[[1]]
    # node_tmp$node_pred = ans[[2]]

    # 判断是否到达了 maximum level，和minimum node size，到达了便退出
    # FACT
    if(sum(node_tmp$portion >= min_nsize) <= 1){
      node_saved[[node_tmp$idx]] = list(node_tmp)
      next # 出口3
    }

    if(node_tmp$misclass > 0){ # If there is no perfect seperation, we choose one of the best seperator
      # 能进这里来的，肯定是可以跑LDA function的
      # We should also record the covariates
      # 不管是什么类型的数据，我们都可以复用

      # Variable selection
      chi_stat = apply(dat_tmp,2,function(x) var_select_all(x,response_tmp, node_tmp$size, Jt, select_method))
      c_split = which(chi_stat == max(chi_stat))[1] # 这个是新的index，不是旧的


      # Splitting
      pjt = prior * node_tmp$portion / Nj
      pjgt = pjt / sum(pjt)
      node_error_rate = sum(pjt) - max(pjt)

      threshold = split_fact(dat_tmp[,c_split],response_tmp,dat_tmp, node_tmp$misclass, pjgt)
      # node_tmp$split_cri = threshold
      # FACT
      if(is.null(threshold)){
        # 用了LDA之后，由于某种原因（大部分是方差过大），导致预测出来全部是同一类
        # 于是决定停止分割
        # node_tmp$split_idx = node_tmp$split_cri = NA
        node_saved[[node_tmp$idx]] = list(node_tmp)
        next # 出口4
      }

      node_tmp$split_cri = round(threshold,1)
      # 下面该判断是否接着分
      no_class = length(node_tmp$split_cri) + 1
      group_idx = cut(dat_tmp[,c_split],breaks = c(-Inf,node_tmp$split_cri,Inf),
                      labels = 1:no_class,right = FALSE) # 先改一下
      idx_children = sapply(1:no_class, FUN = function(o_o) which(group_idx == o_o)) # save the idx_r for children
      subnode_index_c = node_tmp$idx_c # 剩下的cov
      node_tmp$split_idx = node_tmp$idx_c[c_split]

      # 决定是否要继续分了
      children_error_rate = numeric(no_class)
      for(o_o in 1:no_class){
        pjt = prior * table(response_tmp[idx_children[[o_o]]]) / Nj
        # pjgt = pjt / sum(pjt)
        children_error_rate[o_o] = sum(pjt) - max(pjt)
      }

      # 判断是否到达了 maximum level，和minimum node size，到达了便退出
      # FACT
      # cat(node_error_rate,children_error_rate,'\n')
      if(sum(children_error_rate) >= node_error_rate){ # 这里有点问题
        node_tmp$split_idx = node_tmp$split_cri = NA
        node_saved[[node_tmp$idx]] = list(node_tmp)
        next # 出口5
      }

      # If both sides have data (ortherwise, this node should be a leaf node)
      # 下面这个看起来也是一个永远不会遇到的情况

      if(1){
        # Define the spliting rule for numeric variable
        # node_tmp$criteria = ifelse(flag_class, paste(c(c_name[node_tmp$split_idx],'\u2264',
        #                                                ifelse(is.na(node_tmp$split_na_action),'',
        #                                                       ifelse(node_tmp$split_na_action,'**','++')),
        #                                                threshold[1]),collapse = ' '),
        #                            paste(c(c_name[node_tmp$split_idx],'in {', paste(left_group[[1]],collapse = ', '),
        #                                    ifelse(is.na(node_tmp$split_na_action),'',
        #                                           ifelse(node_tmp$split_na_action,'**','++')), '}'),collapse = ' '))
        node_tmp$children = (node_tmp$idx * no_j) : (node_tmp$idx * no_j + o_o -1)
        for(o_o in 1:no_class){
          queue = append(queue,list(list(node_tmp$idx_r[idx_children[[o_o]]],
                                         subnode_index_c, node_tmp$children[o_o])))
        }
      }
    }
    node_saved[[node_tmp$idx]] = list(node_tmp) # 出口4
  }
  # 结果输出部分 ------------------------------------------------------------------

  res = list()
  res$formula = formula
  res$treenode = node_saved # 所有节点信息
  res$dat = dat # 存储源数据——以后要修改，否则占地面积过大
  res$cnames = colnames(dat) # 存储变量名，用来日后的排序
  res$prior = prior
  res$response = response # 为了后面的画图
  cat('The LDA tree is completed.\n')
  return(res)
}
