#' Plot functions
#'
#' Plot the LDA tree (Overall and Single node) using visNetwork and ggplot
#'
#' @param fit A result from LDAtree function
#'
#' @return
#' @export
#' @import visNetwork
#'
#' @examples
plotall <- function(fit){
  ### 全局大图
  node_saved = fit$treenode
  response = fit$response
  stack = c(1)
  idx_curr = 1
  node_count = sum(sapply(node_saved,function(x) !is.null(x))) - 1 # 减去一个 response 的 count
  edges <- data.frame(from = numeric(node_count), to = numeric(node_count))
  id_plot = level_plot = numeric(node_count)
  label_plot = group_plot = character(node_count)

  while(length(stack) > 0){
    id_tmp = node_saved[[stack[1]]][[1]]
    stack = stack[-1]
    id_plot[idx_curr] = id_tmp$idx
    level_plot[idx_curr] = floor(log(id_tmp$idx)/log(2) + 1)
    edges[idx_curr,] = c(floor(id_tmp$idx/2),id_tmp$idx)
    if(!is.na(id_tmp$left)){ # 如果是中间节点，要提供的信息：划分标准，各类占比，总数
      text_zhanbi = paste(round(id_tmp$portion / id_tmp$size,2), collapse = ' / ')
      node_idx_plot = paste('Node',id_tmp$idx)
      label_plot[idx_curr] = paste(id_tmp$criteria, text_zhanbi,
                                   id_tmp$size, node_idx_plot, sep = ' \n ')
      stack = c(stack,id_tmp$left,id_tmp$right)
      group_plot[idx_curr] = NA # 中间节点不分组不涂色
    }else{ # 如果是叶子节点，要提供的信息：预测类别，各类具体，正确/总数
      group_plot[idx_curr] = levels(response)[which.max(id_tmp$portion * fit$prior)]
      text_zhanbi = paste(id_tmp$portion, collapse = ' / ')
      node_idx_plot = paste('Node',id_tmp$idx)
      text_lda = paste(id_tmp$size - id_tmp$misclass, id_tmp$size, sep = ' / ')
      label_plot[idx_curr] = paste(group_plot[idx_curr], text_zhanbi,
                                   text_lda,node_idx_plot, sep = "\n")

    }
    idx_curr = idx_curr + 1L
  }
  edges = edges[-1,]
  nodes <- data.frame(id = id_plot,
                      title = info_panel(fit), # Show when you click
                      level = level_plot,
                      label = label_plot,
                      group = group_plot,
                      shadow = TRUE)

  # 本地预览
  p1 = visNetwork(nodes, edges, width = "100%", height = "400px") %>%
    visNodes(shape = 'dot', color = list(background = "white",
                                         border = "black"))%>%
    visHierarchicalLayout(levelSeparation = 100)%>%
    visLegend(width = 0.1, position = "right", main = "Group")%>%
    visInteraction(dragNodes = FALSE,
                   dragView = TRUE,
                   zoomView = TRUE)
  return(p1)
}

info_panel <- function(fit){
  # 建立信息存储库
  info_node = character(0L)
  node_saved = fit$treenode
  for(i in 1:length(node_saved)){
    node_tmp = node_saved[[i]][[1]]
    if(!is.null(node_tmp)){
      line1 = '#### Information Panel ####'
      line2 = paste('</br>Current Node Index:',node_tmp$idx)
      line3 = paste('</br>There are', node_tmp$size, 'data in this node')
      line4 = paste('</br>The proportion of',paste(levels(fit$response), collapse = ', '),
                    'are', paste(node_tmp$portion, collapse = ', '))
      line5 = paste('</br>There are',node_tmp$misclass, 'error produced by LDA tree')
      error_naive = sum(node_tmp$portion) - max(node_tmp$portion)
      line6 = paste('</br>Compared to stepwise constant tree, the improvement ratio is about',
                    (error_naive - node_tmp$misclass) / error_naive)
      line_final = paste(line1,line2,line3,line4,line5,line6)
      # cat(line_final)
      info_node = c(info_node, line_final)
    }
  }
  return(info_node)
}


