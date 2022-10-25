# GO visualization

# part 1-------
go.vis <- read.csv('../1.part1/20220921_Co down Go pick_By ZQQ.csv',header = T)
go.vis <- go.vis[,-2] 
colnames(go.vis) <- c('Terms','LogP')
go.vis$lp <- -go.vis$LogP
go.vis <- go.vis[order(go.vis$lp),]
go.vis$Terms <- factor(go.vis$Terms,levels = go.vis$Terms)

ggplot(go.vis, aes(x=Terms, y=lp)) +
  geom_segment( aes(x=Terms, xend=Terms, y=0, yend=lp), color="skyblue",size=2,lty = 1) +
  geom_point( color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  xlab("") + ylab(bquote(~-Log[10]~ 'P'))+
  theme(
    axis.text = element_text(size = 12,colour = 'black'),
    axis.title.x = element_text(size = 15,,colour = 'black'),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )
go.vis <- read.csv('../1.part1/20220921_Co up Go pick_By ZQQ.csv',header = T)
go.vis <- go.vis[,-2] 
colnames(go.vis) <- c('Terms','LogP')
go.vis$lp <- -go.vis$LogP
go.vis <- go.vis[order(go.vis$lp),]
go.vis$Terms <- factor(go.vis$Terms,levels = go.vis$Terms)
go.vis <- go.vis[-16,]

ggplot(go.vis, aes(x=Terms, y=lp)) +
  geom_segment( aes(x=Terms, xend=Terms, y=0, yend=lp), color="#e9c46a",size=2) +
  geom_point( color="#e63946", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  xlab("") + ylab(bquote(~-Log[10]~ 'P'))+
  theme(
    axis.text = element_text(size = 12,colour = 'black'),
    axis.title.x = element_text(size = 15,,colour = 'black'),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

# part 4 -----
go.vis <- read.csv('../4.part4/median.cluster15_C14.csv',header = T)
colnames(go.vis) <- c('Terms','LogP')
go.vis$lp <- -go.vis$LogP
go.vis <- go.vis[order(go.vis$lp),]
go.vis$Terms <- factor(go.vis$Terms,levels = go.vis$Terms)

ggplot(go.vis, aes(x=Terms, y=lp)) +
  geom_segment( aes(x=Terms, xend=Terms, y=0, yend=lp), color="#e9c46a",size=2,lty = 1) +
  geom_point( color="#e63946", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  xlab("") + ylab(bquote(~-Log[10]~ 'P'))+
  theme(
    axis.text = element_text(size = 12,colour = 'black'),
    axis.title.x = element_text(size = 15,,colour = 'black'),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

## TF
tf.vis <- read.csv('../4.part4/median.cluster15_C4_14_TF.csv',header = T)
tf.vis$LogP <- -tf.vis$LogP
tf.vis <- tf.vis[order(tf.vis$LogP),]
tf.vis$Description <- factor(tf.vis$Description,levels = tf.vis$Description)

ggplot(data = tf.vis, aes(x = LogP, y = Description)) + 
  geom_point(aes(size = Z.score,color = Enrichment)) +
  scale_color_gradient(low = 'blue',high = 'red3') +
  scale_size(range  =  c(0, 6)) +
  labs( y = 'C4_14 TF targets', x= bquote(~-Log[10]~ 'P')) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10)) 


# part5 --------
go.vis <- read.csv('../5.part5/greenyellow.GOBP.pick.csv',header = T)
colnames(go.vis) <- c('Terms','LogP')
go.vis$lp <- -go.vis$LogP
go.vis <- go.vis[order(go.vis$lp),]
go.vis$Terms <- factor(go.vis$Terms,levels = go.vis$Terms)

ggplot(go.vis, aes(x=Terms, y=lp)) +
  geom_segment( aes(x=Terms, xend=Terms, y=0, yend=lp), color="greenyellow",size=2,lty = 1) +
  geom_point( color="green", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  xlab("") + ylab(bquote(~-Log[10]~ 'P'))+
  theme(
    axis.text = element_text(size = 12,colour = 'black'),
    axis.title.x = element_text(size = 15,,colour = 'black'),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )



