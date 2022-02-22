#' Main function in LDATree
#'
#' Fit the LDATree, and return the treenodes.
#'
#' @param formula an object of class 'formula'
#' @param data Data frame with response variable in it.
#' @param select_method Indicate which methods you want to use to build up the tree.
#' select_method = c('LDATree', 'FACT')
#' @param split_method Indicate how to split the node
#' split_method = c('univariate', 'linear')
#' @param get_size How to get the right size of the tree?
#' get_size = c('CV', 'pre-stopping')
#' @return List
#' @export
#'
#' @examples
LDAtree <- function(){}
Treee <- function(formula, data, select_method = 'FACT', split_method = 'univariate',
                  get_size = 'CV', prior = NULL, max_level = 10, min_nsize = NULL,
                  cv_number = 10, F0 = 4, misclass_cost = NULL){
  # prior 的顺序需要和data中的levels相同
  # data 是含有response的， dat不含有
  message("Good morning, my friend. Wish you a good day :-)")


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
  min_nsize <-  min_nsize %||% max(length(response) %/% 100, 5) # 每个node最少5个
  ###

  ### 设定Prior
  prior <-  prior %||% (tabulate(response)/nrow(data)) # 如果没给prior，用estimated prior
  ###

  ### 设定misclassification cost
  misclass_cost <- misclass_cost %||% (1 - diag(length(levels(response))))
  colnames(misclass_cost) <- colnames(misclass_cost) %||% levels(response)
  rownames(misclass_cost) <- rownames(misclass_cost) %||% levels(response)


  ### 找到所有要使用的X，合在一起变成dat
  dat = data[,sapply(labels(terms(formula, data = data)),function(x) which(x == colnames(data)))]
  ###

  # 找到所有变量的种类
  cov_class = sapply(dat,class) %in% c('numeric', 'integer') # T表示是ordinal

  # FACT 将离散变量转为LD1
  # 现在的处理方式是转成数字之后再也不转回去了
  # 因此我们需要在预测的时候也使用数字
  # 只在画图的时候变成字符
  cat_trans = NULL
  if(pmatch(select_method, c('LDATree', 'FACT')) == 2){
    # cov_class = sapply(dat,class) %in% c('numeric', 'integer')
    if(!all(cov_class)){
      cat_trans = vector("list",dim(dat)[2])
      for(o_o in which(!cov_class)){
        x_x = fact_cat(dat[,o_o],response, prior)
        dat[,o_o] = x_x[[1]] # 拿回变换好的变量
        cat_trans[[o_o]] = x_x[[2]] # 保存变换矩阵
        names(cat_trans)[o_o] = colnames(dat)[o_o]
      }
      cat_trans = cat_trans[!is.na(names(cat_trans))]
    }
  }

  ### 关于NA值的填法 ###
  # 这个得放到FACT变cat为num之后
  # 因为需要拿到的是数字才可以用mle
  NA_method = 'imp_when_needed'
  if(select_method == 'FACT' & split_method == 'linear'){
    NA_method <- 'imp_at_root'
  }

  if(any(sapply(dat,function(o_o) any(is.na(o_o))))){ # 如果存在NA
    if(NA_method == 'imp_at_root'){
      dat <- class_centroid_impute(dat,response, prior, cov_class, cat_trans)
    }
  }
  ###


# 方法预处理 -------------------------------------------------------------------

# 未来这里可以添加很多个参数：
  # 二叉树，还是多叉树


# Model 部分 ----------------------------------------------------------------
  fit = tree_growing(response, dat, prior, misclass_cost, max_level, min_nsize, select_method, split_method, F0, get_size, cov_class, cat_trans)
  # fit = tree_growing_fact(response, dat, prior, max_level, min_nsize, select_method)


  # # FACT 离散变量表
  # if(pmatch(select_method, c('LDATree', 'FACT')) == 2){
  #   if(!is.na(sum(cov_class))){
  #     fit$cat_trans = cat_trans
  #   }
  # }

# Pruning -----------------------------------------------------------------
  if(pmatch(get_size, c('CV', 'pre-stopping')) == 2){
    return(fit) # Pre-stopping rule 直接回到快乐老家
  }
  message('0.0 Now, I am about to start the CV !!!')
  # return(fit)

  if(sum(!sapply(fit$treenode,is.null)) == 1){ # 如果只有根节点，直接出局
    return(fit)
  }

  T_saved = fit$treenode
  fit = traverse(fit) # Get the alpha
  # Divide the data into ten parts
  set.seed(dim(dat)[1]*dim(dat)[2]) # Fixed seed
  idx_CV = sample(c(rep(1:(cv_number-1),each = dim(dat)[1]%/%cv_number),
                    numeric(dim(dat)[1] - (cv_number-1)* (dim(dat)[1]%/%cv_number))+cv_number))
  cv_fit = vector("list", cv_number)
  for(i in 1:cv_number){
    r_tmp = which(idx_CV == i) # 找出第i组的行index
    # 这里有个tree_growing function，需要随时跟着tree_growing的发展来发展
    cv_fit[[i]] = tree_growing(response[-r_tmp], dat[-r_tmp,], prior, misclass_cost, max_level, min_nsize, select_method, split_method, F0, get_size, cov_class, cat_trans) # 找出母树
    ### min_size 好像出了一些问题
    # for(o_o in seq(length(cv_fit[[i]]$treenode))){
    #   print(cv_fit[[i]]$treenode[[o_o]][[1]]$portion)
    # }
    cv_fit[[i]] = traverse(cv_fit[[i]]) # update alpha
  }
  # plotall(cv_fit[[2]])
  CV_table = data.frame(Tree = numeric(0),
                        Tnodes = numeric(0),
                        Mean_MSE = numeric(0),
                        SE_Mean = numeric(0),
                        alpha = numeric(0))
  CV_table[dim(CV_table)[1]+1,] = c(dim(CV_table)[1]+1,sum(!sapply(fit$treenode,is.null)),
                                    get_mean_se(cv_fit,idx_CV,response,dat,cv_number),-1)
  while(sum(!sapply(fit$treenode,is.null)) > 1){ # 除非是切到根结点了
    cat('The current number of node is', sum(!sapply(fit$treenode,is.null)), '\n')
    # 找到切割的alpha
    alpha_tmp = get_cut_alpha(fit)
    fit = cut_alpha(fit,alpha_tmp, select_method)
    fit = traverse(fit)
    for(i in 1:cv_number){
      cv_fit[[i]] = cut_alpha(cv_fit[[i]],alpha_tmp,select_method)
      cv_fit[[i]] = traverse(cv_fit[[i]])
    }
    # plotall(cv_fit[[3]])
    CV_table[dim(CV_table)[1]+1,] = c(dim(CV_table)[1]+1,sum(!sapply(fit$treenode,is.null)),
                                      get_mean_se(cv_fit,idx_CV,response,dat,cv_number),alpha_tmp)
  }
# Evaluation --------------------------------------------------------------
  CV_table = CV_table[dim(CV_table)[1]:1,] # 反向排序，因为which.min会选择最上面的，
  # 但我们需要最小的树
  k_se = 1
  idx_star = which.min(CV_table$Mean_MSE)
  star_threshold = CV_table$Mean_MSE[idx_star] + k_se * CV_table$SE_Mean[idx_star]
  # 2022/02/22 修改：原来的代码不够robust，因为alpha不一定单调
  # 切之后的alpha可能会下降
  # alpha_final = CV_table$alpha[which(CV_table$Mean_MSE<=star_threshold)[1]]
  # 正着数换成反着数，n+1-t
  idx_final = dim(CV_table)[1] + 1 - which(CV_table$Mean_MSE<=star_threshold)[1]
  fit$treenode = T_saved
  fit = traverse(fit) # Get the alpha
  CV_table = CV_table[dim(CV_table)[1]:1,]
  tmp = 1
  # while(CV_table$alpha[tmp] <= alpha_final){
  while(tmp < idx_final){
    fit = cut_alpha(fit,CV_table$alpha[tmp],select_method)
    fit = traverse(fit)
    tmp = tmp + 1
  }
  fit$response_name = colnames(model.frame(formula, data))[1]
  fit$CV_table = CV_table
  cat('The pruned tree is completed\n')
  return(fit)
}



# Tree Growing Function ---------------------------------------------------



tree_growing <- function(response, dat, prior, misclass_cost, max_level, min_nsize, select_method, split_method, F0, get_size, cov_class, cat_trans){
  # Model 部分 ----------------------------------------------------------------
  col_idx = 1:ncol(dat)
  # 对于像FACT来说的多叉树，使用这种方式排列树结构会使得运算效率更高（没人想遍历一个六百万的稀疏列表）
  # idx_next 表示了当前使用的list最后一个位置，也就是长度
  idx_next = 1L
  queue = list(list(1:nrow(dat),col_idx,idx_next,1L, 0L)) # Used to save the current intermediate nodes
  c_name = colnames(dat) # Used to denote the criteria

  node_saved = list() # Save the nodes for printing the tree
  Nj = table(response) # 获得总体proportion，为了计算每个节点的prior
  no_j = length(unique(response)) # 共有多少类


  while(length(queue)!=0){ # 当还有节点没有被划分好

    ### 拿到当前节点的信息
    # if(pmatch(select_method, c('LDATree', 'FACT')) == 1){
    #   # LDATree
    #   node_tmp = generate_node(idx_r = queue[[1]][[1]],
    #                            idx_c = queue[[1]][[2]],
    #                            idx = queue[[1]][[3]]) # Get the current subset index
    # }else if(pmatch(select_method, c('LDATree', 'FACT')) == 2){
    #   # FACT
    #   node_tmp = generate_node_fact(idx_r = queue[[1]][[1]],
    #                                 idx_c = queue[[1]][[2]],
    #                                 idx = queue[[1]][[3]]) # Get the current subset index
    # }
    node_tmp = generate_node(idx_r = queue[[1]][[1]],
                               idx_c = queue[[1]][[2]],
                               idx = queue[[1]][[3]],
                             layer = queue[[1]][[4]],
                             parent = queue[[1]][[5]], select_method) # Get the current subset index

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
    node_tmp$pred_method = 'mode' # 默认用众数来估计
    # node_tmp$node_pred = levels(response)[which.max(node_tmp$portion * prior)] # 这个之后反正可以再改
    # node_tmp$node_pred = levels(response)[which.min(misclass_cost %*% (node_tmp$portion * prior))] # 这个之后反正可以再改
    node_tmp$node_pred = mis_cost_cal(proportion = node_tmp$portion * prior, misclass_cost = misclass_cost,
                                     method = 'pred', level = levels(response_tmp)) # 这个之后反正可以再改

    ###

    ### 出口1: 无X，或者所有Y都相同
    Jt = sum(node_tmp$portion != 0) # number of class in node t 当前节点有几类
    if(Jt == 1 | node_tmp$covs == 0){ # If all of the data belongs to one class, or there is no covariates left
      node_tmp$misclass = node_tmp$size - max(node_tmp$portion) # 既然分不下去了，那就只能预测众数
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
      node_saved[[node_tmp$idx]] = list(node_tmp)
      next # 出口2
    }

    # fit <- total_LDA(dat_tmp, response_tmp)
    # predict_tmp = predict(fit,newdata = dat_tmp) # Predict using the fitted model
    # node_tmp$misclass = sum(predict_tmp$class != response_tmp)

    # node_tmp$misclass = get_error_LDA(dat_tmp, response_tmp, prior) # 这个信息会一直保留，因为最后需要展示在图上
    if(pmatch(select_method, c('LDATree', 'FACT')) == 1){
      # LDATree需要在节点处fitLDA以进行错误计算
      ans = pred_LDA(dat_tmp, response_tmp, prior)
      node_tmp$misclass = ans[[3]] # 这个信息会一直保留，因为最后需要展示在图上
      node_tmp$pred_method = ans[[1]]
      node_tmp$node_pred = ans[[2]]
    }else if(pmatch(select_method, c('LDATree', 'FACT')) == 2){
      # FACT
      # node_tmp$misclass = node_tmp$size - max(node_tmp$portion)
      # 这个misclass暂时先不加prior了
      node_tmp$misclass <- mis_cost_cal(proportion = node_tmp$portion, misclass_cost = misclass_cost,
                   level = levels(response_tmp))
    }


    # 判断是否到达了 maximum level，和minimum node size，到达了便退出
    if(pmatch(select_method, c('LDATree', 'FACT')) == 1){
      # LDATree
      if(node_tmp$idx >= 2 ** (max_level-1) | node_tmp$size < min_nsize){
        # ans = pred_LDA(dat_tmp, response_tmp, prior)
        # node_tmp$pred_method = ans[[1]]
        # node_tmp$node_pred = ans[[2]]
        node_saved[[node_tmp$idx]] = list(node_tmp)
        next # 出口3
      }
    }else if(pmatch(select_method, c('LDATree', 'FACT')) == 2){
      # FACT
      # 如果最多只有一组大于最小的node size
      # cat('The current layer is:', node_tmp$layer,'\n')
      if(node_tmp$layer >= max_level | (sum(node_tmp$portion >= min_nsize) <= 1)){
        node_saved[[node_tmp$idx]] = list(node_tmp)
        next # 出口3
      }
    }




    if(node_tmp$misclass > 0){ # If there is no perfect seperation, we choose one of the best seperator
      # 能进这里来的，肯定是可以跑LDA function的
      # We should also record the covariates
      # 不管是什么类型的数据，我们都可以复用


# Variable selection ------------------------------------------------------



      if(pmatch(select_method, c('LDATree', 'FACT')) == 2){
        # FACT
        if(pmatch(split_method, c('univariate', 'linear')) == 1){
          # 单变量的变量选择
          # FACT

          # Variable Selection
          F_stat = apply(dat_tmp,2,function(x) var_select_all(x,response_tmp, node_tmp$size, Jt, select_method))
          c_split = which.max(F_stat) # 这个是新的index，不是旧的



          if(F_stat[c_split] < F0){
            # Change it to dispersion
            dat_tmp_trans = abs(scale(dat_tmp,center = TRUE, scale = FALSE))
            F_stat_trans = apply(dat_tmp_trans,2,function(x) var_select_all(x,response_tmp, node_tmp$size, Jt, select_method))
            c_split_trans = which.max(F_stat_trans)
            x_split_mean = attr(dat_tmp_trans,'scaled:center')[c_split_trans] # 保存一下均值
            if(F_stat_trans[c_split_trans] >= F0){
              # 如果Zi显著
              node_tmp$dispersion = 1 # 记录下来mean
              x_split_save = dat_tmp[,c_split_trans]
              dat_tmp[,c_split_trans] = dat_tmp_trans[,c_split_trans]
              F_stat = F_stat_trans
              c_split = c_split_trans
            }else{
              # 用原本的mean来划分
              node_tmp$dispersion = 2
            }
          }



          # prior calculation
          pjt = prior * node_tmp$portion / Nj
          pjgt = pjt / sum(pjt) # standardize
          node_error_rate = mis_cost_cal(proportion = pjt, misclass_cost = misclass_cost) # 这里要用绝对错误，而不是相对错误



          ### 检查是否具有缺失值 ###
          if(any(is.na(dat_tmp[,c_split]))){
            # 这里对于categorical变量我们没有进行离散化处理
            # 完全是出于写代码的便利，加上我认为是否离散化区别也不大
            # 感兴趣的朋友可以自己改
            dat_tmp <- class_centroid_impute(dat_tmp, y, prior, cov_class, cat_trans)
          }

          threshold = split_fact_uni(dat_tmp[,c_split],response_tmp, pjgt, misclass_cost, min_nsize)
          # node_tmp$split_cri = threshold
          # FACT
          if(is.null(threshold)){
            # 用了LDA之后，由于某种原因（大部分是方差过大），导致预测出来全部是同一类
            # 于是决定停止分割
            # node_tmp$split_idx = node_tmp$split_cri = NA
            node_saved[[node_tmp$idx]] = list(node_tmp)
            next # 出口4
          }

          node_tmp$split_cri = round(threshold,2) # 只保留一位小数
          # node_tmp$split_cri = threshold # 只保留一位小数

          # 下面该判断是否接着分
          # 如果是Dispersion来的，我们需要将原来的数据变回去，同时修改split
          # 这下面一通操作，我已经看不懂了
          if(!is.na(node_tmp$dispersion)){
            if(node_tmp$dispersion == 2L){
              node_tmp$split_cri = x_split_mean
              cat('The dispersion with mean is runing!\n')
            }else{
              cat('The dispersion with multi-split is runing!\n')
              dat_tmp[,c_split] = x_split_save
              final_cut = sort(unique(c(x_split_mean-node_tmp$split_cri, x_split_mean + node_tmp$split_cri)))

              # 调试的时候出现了一个很神奇的bug
              # 本来是满足最小节点数的，但是因为dispersion会因为绝对值的关系
              # 把一个分割拆成两个，就导致不满足最小节点数，所以这里加上这个判别条件
              # 下面的代码复制于split_fact_uni
              final_cut_pro = final_cut
              o_o = 1
              stop_split_flag = 0
              while(o_o <= length(final_cut_pro)){
                if(o_o == 1){
                  # 初始化分组信息
                  test_tmp = as.numeric(cut(dat_tmp[,c_split],breaks = c(-Inf,final_cut,Inf))) # 强制变成integer
                }
                # 如果最大组不大于最低要求的话，要和右面的组合并
                # print(table(y[which(test_tmp == o_o)]))
                if(sum(table(response_tmp[which(test_tmp == o_o)]) >= min_nsize) < 1){
                  if(length(final_cut_pro) == 1){ # 如果被切的一个也不剩了，停止划分
                    stop_split_flag = 1
                    break
                  }
                  final_cut_pro = final_cut_pro[-o_o]
                  o_o = 1
                }else{
                  o_o = o_o + 1
                }
              }

              # 实践中发现最右面的组有可能不满足min_size的要求
              if(sum(table(response_tmp[which(test_tmp == o_o)]) >= min_nsize) < 1){
                if(length(final_cut_pro) == 1){ # 如果被切的一个也不剩了，停止划分
                  stop_split_flag = 1
                }
                final_cut_pro = final_cut_pro[-length(final_cut_pro)]
              }

              if(stop_split_flag){
                node_saved[[node_tmp$idx]] = list(node_tmp)
                next # 出口4.1
              }
              node_tmp$split_cri = final_cut_pro

              # final_cut_pro = c()
              #
              # test_tmp = table(cut(dat_tmp[,c_split],breaks = c(-Inf,final_cut,Inf)))
              # for(o_o in 1 : length(final_cut)){
              #   if(test_tmp[o_o] != 0){
              #     final_cut_pro = c(final_cut_pro, final_cut[o_o])
              #   }
              # }
              # if(test_tmp[length(test_tmp)] == 0){
              #   final_cut_pro = final_cut_pro[-length(final_cut_pro)]
              # }
              # node_tmp$split_cri = final_cut_pro
            }
          }

          no_class = length(node_tmp$split_cri) + 1 #  产生了几个子节点
          group_idx = cut(dat_tmp[,c_split],breaks = c(-Inf,node_tmp$split_cri,Inf),
                          labels = 1:no_class,right = TRUE) # 包含右边界
          idx_children = sapply(1:no_class, FUN = function(o_o) which(group_idx == o_o)) # save the idx_r for children
          subnode_index_c = node_tmp$idx_c # 剩下的cov
          node_tmp$split_idx = node_tmp$idx_c[c_split]

          # 这段代码要保留
          # 判断是否不存在任何一组大于minimum node size，到达了便退出
          # ephemeral = sapply(seq(no_class), function(o_o) sum(table(response_tmp[idx_children[[o_o]]]) >= min_nsize))
          # if(min(ephemeral) < 1){
          #   # print(table(response_tmp[idx_children[[o_o]]]))
          #   # cat('NA happens there 1:', node_tmp$idx,'\n')
          #   print('Exit X')
          #   node_tmp$split_idx = node_tmp$split_cri = NA
          #   node_saved[[node_tmp$idx]] = list(node_tmp)
          #   next # 出口X
          # }



          # 决定是否要继续分了
          # 只有pre-stopping会根据错误率是否下降来判断停止与否
          if(pmatch(get_size, c('CV', 'pre-stopping')) == 2){
            children_error_rate = numeric(no_class)
            for(o_o in 1:no_class){
              pjt = prior * table(response_tmp[idx_children[[o_o]]]) / Nj
              # pjgt = pjt / sum(pjt)
              children_error_rate[o_o] = mis_cost_cal(proportion = pjt, misclass_cost = misclass_cost)
            }
            if(sum(children_error_rate) >= node_error_rate){ # 这里有点问题 # ？？？啥问题
              # cat('NA happens there 2:', node_tmp$idx,'\n')
              node_tmp$split_idx = node_tmp$split_cri = NA
              node_saved[[node_tmp$idx]] = list(node_tmp)
              next # 出口5
            }
          }
          ##### 这里是单变量划分的结束 #####
        }else if(pmatch(split_method, c('univariate', 'linear')) == 2){
          # 通过PCA拿到一堆Y
          cat('The dimension of the dat_tmp:', dim(dat_tmp)[2],'\n')
          res_tmp = PCA_Y(dat_tmp, response_tmp)
          dat_tmp_y = res_tmp[[1L]] # 这时候剩几维我们已经不好讲了
          linear_split_trans <- res_tmp[[2L]] # for prediction
          #对这些Y做变量选择
          F_stat = apply(dat_tmp_y,2,function(x) var_select_all(x,response_tmp, node_tmp$size, Jt, select_method))
          c_split = which.max(F_stat)
          # if(F_stat[c_split] >= F0){
          if(FALSE){
            # 进行正常的纯linear split
            # Situation: 1
            cat('Linear Combination Scenario 1\n')
            ephemeral = fact_univar(dat_tmp_y[,c_split], response_tmp, node_tmp, prior, misclass_cost, Nj, min_nsize, get_size)
            if(is.null(ephemeral)){
              node_saved[[node_tmp$idx]] = list(node_tmp)
              next # 出口
            }else{
              # 设置划分
              node_tmp = ephemeral[[1]]
              node_tmp$split_idx = node_tmp$idx_c[c_split]
              node_tmp$linear_split_trans = linear_split_trans
              # 这个地方肯定是跟environment有点关系，我发现如果写在外面包成函数的话，
              # 所有节点的linear split function 都一模一样，就像我写了x_new[,node_tmp$idx_c]
              # 之后，所有节点最后都取的是最后一个节点的idx_c，很怪
              idx_children = ephemeral[[2]]
              no_class = ephemeral[[3]]
            }
            ####
          }else{
            trans_w_tmp = scale(dat_tmp_y,center = TRUE, scale = FALSE)
            dat_tmp_w = abs(trans_w_tmp)
            F_stat2 = apply(dat_tmp_w,2,function(x) var_select_all(x,response_tmp, node_tmp$size, Jt, select_method))
            c_split2 = which.max(F_stat2)
            if(F_stat2[c_split2] < F0){
              # univariate split, on y average
              # Situation: 2
              cat('Linear Combination Scenario 2\n')
              ephemeral = fact_univar(dat_tmp_y[,c_split], response_tmp, node_tmp, prior, misclass_cost, Nj, min_nsize, get_size, simple_mean = TRUE)
              if(is.null(ephemeral)){
                node_saved[[node_tmp$idx]] = list(node_tmp)
                next # 出口
              }else{
                node_tmp = ephemeral[[1]]
                node_tmp$linear_split_trans = linear_split_trans
                node_tmp$split_idx = node_tmp$idx_c[c_split]
                idx_children = ephemeral[[2]]
                no_class = ephemeral[[3]]
              }

            }else{
              # levene's test on x
              p_stat = p.adjust(apply(dat_tmp_w,2,function(o_o) car::leveneTest(y = o_o, group = response_tmp)[[3]][1]))
              levene_sig <- (p_stat <= 0.05)
              if(sum(levene_sig) == 1){
                # univariate split on w
                # Situation: 3
                cat('Linear Combination Scenario 3\n')
                # print(levene_sig)
                # print(dim(dat_tmp_y))
                # print(dim(dat_tmp_w))
                ephemeral = fact_univar(dat_tmp_w[,levene_sig], response_tmp, node_tmp, prior, misclass_cost, Nj, min_nsize, get_size)
                if(is.null(ephemeral)){
                  node_saved[[node_tmp$idx]] = list(node_tmp)
                  next # 出口
                }else{
                  node_tmp = ephemeral[[1]]
                  node_tmp$split_idx = node_tmp$idx_c[c_split]
                  node_tmp$linear_split_trans = PCA_Y(dat_tmp, response_tmp, situation = 3)
                  idx_children = ephemeral[[2]]
                  no_class = ephemeral[[3]]
                }

              }else if(sum(levene_sig) == 0){
                # change everything to polar
                # univariate split on r,theta
                # Situation: 4
                cat('Linear Combination Scenario 4\n')
                dat_tmp_polar = nsphere(dat_tmp_w)
                F_stat = apply(dat_tmp_polar,2,function(x) var_select_all(x,response_tmp, node_tmp$size, Jt, select_method))
                c_split = which.max(F_stat)
                ephemeral = fact_univar(dat_tmp_polar[,c_split], response_tmp, node_tmp, prior, misclass_cost, Nj, min_nsize, get_size)
                if(is.null(ephemeral)){
                  node_saved[[node_tmp$idx]] = list(node_tmp)
                  next # 出口
                }else{
                  node_tmp = ephemeral[[1]]
                  node_tmp$split_idx = node_tmp$idx_c[c_split]
                  node_tmp$linear_split_trans = PCA_Y(dat_tmp, response_tmp, situation = 4)
                  idx_children = ephemeral[[2]]
                  no_class = ephemeral[[3]]
                }

              }else{
                # change only the significant ones
                # univariate split on w,r,theta
                # Situation: 5
                cat('Linear Combination Scenario 5\n')
                dat_tmp_polar = dat_tmp_w
                dat_tmp_polar[,levene_sig] = nsphere(dat_tmp_polar[,levene_sig])
                F_stat = apply(dat_tmp_polar,2,function(x) var_select_all(x,response_tmp, node_tmp$size, Jt, select_method))
                c_split = which.max(F_stat)
                ephemeral = fact_univar(dat_tmp_polar[,c_split], response_tmp, node_tmp, prior, misclass_cost, Nj, min_nsize, get_size)
                if(is.null(ephemeral)){
                  node_saved[[node_tmp$idx]] = list(node_tmp)
                  next # 出口
                }else{
                  node_tmp = ephemeral[[1]]
                  node_tmp$split_idx = node_tmp$idx_c[c_split]
                  node_tmp$linear_split_trans = PCA_Y(dat_tmp, response_tmp, situation = 5)
                  idx_children = ephemeral[[2]]
                  no_class = ephemeral[[3]]
                }
              }
            }
          }
        }
      }



# Variable Selection & Splitting ---------------------------------------------------------------



      if(pmatch(select_method, c('LDATree', 'FACT')) == 1){
        # LDATree

        # Variable Selection
        chi_stat = apply(dat_tmp,2,function(x) var_select_all(x,response_tmp, node_tmp$size, Jt, select_method))
        c_split = which.max(chi_stat) # 这个是新的index，不是旧的

        flag_class = class(dat_tmp[,c_split]) %in% c('numeric', 'integer')
        if(flag_class){
          threshold = split_noncat(dat_tmp[,c_split],response_tmp,dat_tmp, node_tmp$misclass, prior)
          node_tmp$split_cri = threshold[1]
          node_tmp$split_na_action = threshold[2]
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
        ####### 结束 #######
      }


# Generating Children Nodes -----------------------------------------------

      if(pmatch(select_method, c('LDATree', 'FACT')) == 1){
        # LDATree
        if(length(idx_left) != 0 & length(idx_right) != 0){
          # Define the spliting rule for numeric variable
          node_tmp$criteria = ifelse(flag_class, paste(c(c_name[node_tmp$split_idx],'\u2264',
                                                         ifelse(is.na(node_tmp$split_na_action),'',
                                                                ifelse(node_tmp$split_na_action,'**','++')),
                                                         threshold[1]),collapse = ' '),
                                     paste(c(c_name[node_tmp$split_idx],'in {', paste(left_group[[1]],collapse = ', '),
                                             ifelse(is.na(node_tmp$split_na_action),'',
                                                    ifelse(node_tmp$split_na_action,'**','++')), '}'),collapse = ' '))

          node_tmp$left = 2 * node_tmp$idx # Decide the child index for the current node
          queue = append(queue,list(list(node_tmp$idx_r[idx_left],subnode_index_c, node_tmp$left, node_tmp$layer + 1, node_tmp$idx)))
          node_tmp$right = 2 * node_tmp$idx + 1 # Decide the child index for the current node
          queue = append(queue,list(list(node_tmp$idx_r[idx_right],subnode_index_c, node_tmp$right, node_tmp$layer + 1, node_tmp$idx)))
          # 为了pruning，要把left and right 都放进 children
          node_tmp$children = c(node_tmp$left, node_tmp$right)

        }
      }else if(pmatch(select_method, c('LDATree', 'FACT')) == 2){
        # FACT
        # cat('o_o is:', o_o, '\n')
        # node_tmp$children = (node_tmp$idx * 2 * no_j) : (node_tmp$idx * 2 * no_j + no_class -1)
        # print(node_tmp$children)
        node_tmp$children = seq(no_class) + idx_next
        for(o_o in 1:no_class){
          queue = append(queue,list(list(node_tmp$idx_r[idx_children[[o_o]]],
                                         node_tmp$idx_c, node_tmp$children[o_o], node_tmp$layer + 1, node_tmp$idx)))
        }
        idx_next = idx_next + no_class # 刷新一下maximum队列
        # cat('idx_next is:', idx_next, '\n')
      }
    }
    # cat('The index is:',node_tmp$idx,'\n')
    # print(list(node_tmp))
    node_saved[[node_tmp$idx]] = list(node_tmp) # 出口4
  }
  # 结果输出部分 ------------------------------------------------------------------

  res = list()
  res$formula = formula
  res$select_method = select_method
  res$treenode = node_saved # 所有节点信息
  # linear split的话，所有位置都是数字
  res$cat_trans = cat_trans
  # res$dat = dat # 存储源数据——以后要修改，否则占地面积过大
  res$cnames = colnames(dat) # 存储变量名，用来日后的排序
  res$prior = prior
  res$response = response # 为了后面的画图
  # prepare for the imformation in output
  res$split_method = split_method # for prediction
  res$misclass_cost = misclass_cost

  # 2022/02/22 下面这个代码看得我一头雾水，what is this?
  # My guess: somehow univariate 需要全是TRUE才能运行
  if(select_method == 'FACT' & split_method == 'univariate'){
    res$cov_class = cov_class > -2
  }else{
    res$cov_class = cov_class
  }

  cat('The LDA tree is completed.\n')
  return(res)
}


