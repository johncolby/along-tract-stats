exDir    = '/path/to/along-tract-stats/example'
subID    = 'subject1'
trk_info = read.table(file.path(exDir, 'tract_info.txt'), header=T, sep='\t')
trk_avgs = read.table(file.path(exDir, 'trk_avgs.txt'), header=T, sep='\t')

library(ggplot2)

################################################################################
# Loop over tracts
for(iTrk in 1:nrow(trk_info)){
    # Import/format data
    trkName = sprintf('%s_%s', trk_info$Tract[iTrk], trk_info$Hemisphere[iTrk])
    path = file.path(exDir, sprintf('%s_%s.txt', subID, trkName))
    single_sub            = read.table(path, header=T)
    single_sub$Streamline = factor(single_sub$Streamline)
    single_sub$Point      = factor(single_sub$Point)
    
    # Fit linear model
	cat(paste('Tract:', as.character(trk_info$Tract[iTrk])))
	cat('\n')
	print(anova(lm(FA ~ Point, data=single_sub)))
	cat('\n')
    
    # Plot each streamline (FA vs. position)
    # Address overplotting with alpha and a slight x-position jitter
    scale     = 100/(max(as.numeric(single_sub$Point))-1)
    p = ggplot(single_sub, aes(x=(as.numeric(Point)-1)*scale, y=FA, group=Streamline)) + geom_line(alpha=0.01, position='jitter', width=0.01)
    
    # Overlay along-tract mean FA ± SD
    get_band <- function(y.in){
        ymax = mean(y.in)+sd(y.in)
        ymin = mean(y.in)-sd(y.in)
        data.frame(ymax,ymin)
    }
    means = c(stat_summary(aes(group=1), fun.y=mean, geom='line', color='blue'), stat_summary(aes(group=1), geom='ribbon', fun.data=get_band, color=0, fill='blue', alpha=0.2))
    
    # Overlay tract-averaged mean FA ± SD
    fa.avg = trk_avgs$Mean.FA[iTrk]
    fa.sd  = trk_avgs$SD.FA[iTrk]
    avg = geom_pointrange(aes(x=50, y=fa.avg, ymax=fa.avg+fa.sd, ymin=fa.avg-fa.sd), color='red')
    
    # Options
    opts = list(xlab('Position along tract (%)'), coord_cartesian(xlim = (c(-5,105))), ylim(c(0,1)))
    
    # Draw final plot
    dev.new(width=4, height=4)
    print(p + means + avg + opts)
    #ggsave(file.path(exDir, sprintf('%s.jpg', trkName)), dpi=600)
}