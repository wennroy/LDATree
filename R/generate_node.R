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
generate_node <- function(idx_r = NA,idx_c = NA, idx = NA){ #还剩哪些行，哪些列？
  # 类的创建
  me <- list(
    idx = idx, # 第几个节点
    idx_r = idx_r, # Row
    idx_c = idx_c, # Column
    size = length(idx_r), # 有多少行
    covs = length(idx_c), # 有多少列
    left = NA, # 左孩子节点
    right = NA, # 右孩子节点
    misclass = NA, # 有多少个分类错误
    portion = NA, # 每一类有多少个数据
    # alpha = NA, # for CART pruning
    leaves = c(), # 所有后代
    criteria = NA, # 用来打印在output tree上面
    split_idx = NA, # 用来记录是用哪一个变量进行split
    split_cri = NA, # Splitting criteria
    pred_method = NA,
    # prior = NA, # group prior for each class, 目前感觉没啥用
    lda_pred = NA # Function 用来预测新进来的数据如果不走了，原地预测
  )
  # Set the name for the class
  class(me) <- append(class(me), "Treenode")
  return(me)
}
