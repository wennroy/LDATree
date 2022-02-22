#' Tree node structure
#'
#' @param idx_r
#' @param idx_c
#' @param idx
#'
#' @return
#' @export
#'
#' @examples
generate_node <- function(idx_r = NA,idx_c = NA, idx = NA, layer = NA, parent = NA, select_method){ #还剩哪些行，哪些列？
  # 类的创建
  me <- list(
    idx = idx, # 第几个节点
    idx_r = idx_r, # Row
    idx_c = idx_c, # Column
    size = length(idx_r), # 有多少行
    covs = length(idx_c), # 有多少列
    # left = NA, # 左孩子节点
    # right = NA, # 右孩子节点
    misclass = NA, # 有多少个分类错误
    portion = NA, # 每一类有多少个数据
    alpha = NA, # for CART pruning
    leaves = c(), # 所有后代，可以用是否是NULL来判断是否为叶子结点
    layer = layer, # 用来记录当前的层数
    parent = parent, # 用来记录父亲是谁
    children = NA, # Children 本来是多叉树的特征
    criteria = NA, # 用来打印在output tree上面
    split_idx = NA, # 用来记录是用哪一个变量进行split，也可以用来判断是否为叶子结点
    split_cri = NA, # Splitting criteria
    split_na_action = NA, # 把NA分到左面或者右面, 1 代表左面，0代表右面。NA代表训练时没有NA
    pred_method = NA, #
    linear_split_trans = NA, # for linear combination split, record the reverse function for prediction
    # prior = NA, # group prior for each class, 目前感觉没啥用
    resub_error = NA, # R(T): for CV pruning
    node_pred = NA # Function 用来预测新进来的数据如果不走了，原地预测
  )
  # 孩子节点的个数：几叉树
  if(pmatch(select_method, c('LDATree', 'FACT')) == 1){
    # LDATree
    me$left = NA # 左孩子节点
    me$right = NA # 右孩子节点
  }else if(pmatch(select_method, c('LDATree', 'FACT')) == 2){
    # FACT
    me$dispersion = NA
  }
  # Set the name for the class
  class(me) <- append(class(me), "Treenode")
  return(me)
}

# generate_node_fact <- function(idx_r = NA,idx_c = NA, idx = NA){ #还剩哪些行，哪些列？
#   # 类的创建
#   me <- list(
#     idx = idx, # 第几个节点
#     idx_r = idx_r, # Row
#     idx_c = idx_c, # Column
#     size = length(idx_r), # 有多少行
#     covs = length(idx_c), # 有多少列
#     children = NA, # 孩子节点 FACT
#     misclass = NA, # 有多少个分类错误
#     portion = NA, # 每一类有多少个数据
#     alpha = NA, # for CART pruning
#     leaves = c(), # 所有后代，可以用是否是NULL来判断是否为叶子结点
#     criteria = NA, # 用来打印在output tree上面
#     split_idx = NA, # 用来记录是用哪一个变量进行split，也可以用来判断是否为叶子结点
#     split_cri = NA, # Splitting criteria
#     split_na_action = NA, # 把NA分到左面或者右面, 1 代表左面，0代表右面。NA代表训练时没有NA
#     pred_method = NA, #
#     # prior = NA, # group prior for each class, 目前感觉没啥用
#     resub_error = NA, # R(T): for CV pruning
#     node_pred = NA # Function 用来预测新进来的数据如果不走了，原地预测
#   )
#   # Set the name for the class
#   class(me) <- append(class(me), "Treenode")
#   return(me)
# }
