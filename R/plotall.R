#' Plot functions
#'
#' Plot the LDA tree (Overall and Single node) using visNetwork and ggplot
#'
#' @param fit A result from LDAtree function
#'
#' @return
#' @export
#'
#' @examples
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

# LDA 可视化 -----------------------------------------------------------------


### 全局大图

plotall <- function(fit){
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
      label_plot[idx_curr] = paste(id_tmp$criteria, text_zhanbi, id_tmp$size, sep = ' \n ')
      stack = c(stack,id_tmp$left,id_tmp$right)
      group_plot[idx_curr] = NA # 中间节点不分组不涂色
    }else{ # 如果是叶子节点，要提供的信息：预测类别，各类具体，正确/总数
      group_plot[idx_curr] = levels(response)[which.max(id_tmp$portion)]
      text_zhanbi = paste(id_tmp$portion, collapse = ' / ')
      text_lda = paste(id_tmp$size - id_tmp$misclass, id_tmp$size, sep = ' / ')
      label_plot[idx_curr] = paste(group_plot[idx_curr], text_zhanbi, text_lda, sep = "\n")

    }
    idx_curr = idx_curr + 1L
  }
  edges = edges[-1,]
  nodes <- data.frame(id = id_plot,
                      title = info_panel(fit), # Show when you click
                      level = level_plot,
                      label = label_plot,
                      group = group_plot)

  # 本地预览
  p1 = visNetwork::visNetwork(nodes, edges, width = "100%", height = "400px") %>%
    visNetwork::visNodes(shape = 'dot', color = list(background = "white",
                                         border = "black"))%>%
    visNetwork::visHierarchicalLayout(levelSeparation = 100)%>%
    visNetwork::visLegend(width = 0.1, position = "right", main = "Group")%>%
    visNetwork::visInteraction(dragNodes = FALSE,
                   dragView = TRUE,
                   zoomView = TRUE)
  return(p1)
}

### 局部小图

plotsingle <- function(fit, idx){
  node_tmp = fit$treenode[[idx]][[1]]
  idx_r = node_tmp$idx_r
  idx_c = node_tmp$idx_c
  response_tmp = fit$response[idx_r]

  # 当数据中只有一种response的时候，画出一个竖着的图
  if(length(unique(response_tmp)) == 1){
    dat_new = data.frame(x = numeric(node_tmp$size), y = runif(node_tmp$size),
                         response = response_tmp, pred = numeric(node_tmp$size) + 1.5)
    p = ggplot(data = dat_new)+
      geom_point(aes(x = x, y = y, color = response),size = dat_new$pred, alpha = 0.5)+
      labs(x = 'LD1',y = 'LD2')+
      theme_bw()+
      scale_color_discrete(name=fit$response_name)
    # theme(text=element_text(size=23))
    return(p)
  }

  dat_combined = cbind(fit$dat[idx_r, idx_c],response_tmp)
  colnames(dat_combined)[node_tmp$covs+1] = fit$response_name

  idx_same = apply(dat_combined,2,function(x) length(unique(x)) == 1)
  if(any(idx_same)){
    dat_combined = dat_combined[,-which(idx_same)] # 删除那些x只有一个值的变量
  }
  fit_new = MASS::lda(fit$formula, data = dat_combined) # fit好了的LDA模型

  # 将LD1 & LD2保留，顺便将Y变成predicted，这样方便contour画出图
  # 一本书的作者说， 边界不是理论出来的而是带数试出来的，所以可能会弯弯曲曲

  # 如果只有1个LD1
  LD_df = predict(fit_new,newdata = dat_combined)$x
  # ddd <- data.frame(LD1 = LD_df[,1], LD2 = runif(nrow(LD_df)), group = response_tmp)
  # ggplot(data = ddd)+
  #   geom_point(aes(x = LD1, y = LD2, pch = response_tmp, color = response_tmp), size = 3)

  if(ncol(LD_df) == 1){
    set.seed(round(log(nrow(LD_df))) * 2) # 选一个不变的常数
    LD_df = data.frame(LD_df, LD2 = runif(nrow(LD_df)))
  }
  # datPred <- data.frame(pred = predict(fit_new)$class,LD_df[,1:2])
  datPred <- data.frame(pred = response_tmp, LD_df[,1:2])

  #Create decision boundaries
  fit2 <- MASS::lda(pred ~ LD1 + LD2, data= datPred)

  ld1lim <- scales::expand_range(range(datPred$LD1),mul=0.05)
  ld2lim <- scales::expand_range(range(datPred$LD2),mul=0.05)
  ld1 <- seq(ld1lim[1], ld1lim[2], length.out=300)
  ld2 <- seq(ld2lim[1], ld2lim[2], length.out=300)

  newdat <- expand.grid(list(LD1=ld1,LD2=ld2)) # 拼出了一个300*300的网格数据
  preds <-predict(fit2,newdata=newdat)

  df <- data.frame(x=newdat$LD1, y=newdat$LD2, class = preds$class)
  df$classnum <- as.numeric(df$class)

  # 这个函数很强，输入数字，就可以告诉你如果是前K组，那么ggplot会用什么颜色
  # colorfun <- function(n,l=65,c=100){
  #   hues = seq(15, 375, length=n+1)
  #   return(hcl(h=hues, l=l, c=c)[1:n])
  # } # default ggplot2 colours
  # colors <- colorfun(4)
  p = ggplot(datPred, aes(x=LD1, y=LD2)) +
    geom_raster(data=df, aes(x=x, y=y, fill = factor(class)),alpha=0.2, show.legend=FALSE)+ # 背景填色
    geom_contour(data=df, aes(x=x, y=y, z=classnum), colour="black", alpha=0.5, size = 0.1)+ # 划出分割线
    geom_point(data = datPred, size = 3, aes(pch = pred,colour=pred), alpha = 0.7)+
    theme_bw()+
    scale_x_continuous(expand=c(0,0))+ # 为了使x,y的边界不会有一段突兀的空白
    scale_y_continuous(expand=c(0,0))+
    scale_color_discrete(name=fit$response_name) +
    scale_shape_discrete(name=fit$response_name)
  return(p)
}


