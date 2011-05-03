exptDir = '/path/to/along-tract-stats/example/between_group'
grpLabs = c('Control', 'FASD')
thresh  = 0.05
nPerms  = 100

library(nlme)         # Mixed-effects models
library(ggplot2)      # Plotting tools
library(plyr)         # Data manipulation
library(RColorBrewer) # Color tables

################################################################################
# Import and format data
# Read in demographics
demog       = read.table(file.path(exptDir, 'Demographics.txt'), header=T)
demog$Group = factor(demog$Group, levels=rev(levels(demog$Group)), labels=grpLabs)

# Read in whole-track properties (ex: streamlines) and merge with demographics
trk_props_long = read.table(file.path(exptDir, 'trk_props_long.txt'), header=T)
trk_props_long = merge(trk_props_long, demog)

# Read in length-parameterized track data (ex: FA) and merge with demographics
trk_data              = read.table(file.path(exptDir, 'trk_data.txt'),  header=T)
trk_data$Point        = factor(trk_data$Point)
trk_data[trk_data==0] = NA
trk_data              = merge(trk_data, trk_props_long)
# Add a Position column for easier plotting
trk_data              = ddply(trk_data, c("Tract", "Hemisphere"),
                        transform, Position = (as.numeric(Point)-1) * 100/(max(as.numeric(Point))-1))

################################################################################
# Fit LME models
# Overall ANOVA for Group and Point:Group effects
fit_trk_model1 <- function(df){
    lme.trk = lme(FA ~ Point*Group, data=df, random = ~ 1 | ID, na.action=na.omit)
    data.frame(Term = rownames(anova(lme.trk)), anova(lme.trk))
}
models = list()
models$anova = ddply(trk_data, c("Tract", "Hemisphere"), fit_trk_model1)

# Fit a cell-means version to get effect sizes relative to controls
fit_trk_model2 <- function(df){
    lme.trk = tryCatch(lme(FA ~ Point/Group - 1, data=df, random = ~ 1 | ID, na.action=na.omit), error = function(e) data.frame())
    if(length(lme.trk)!=0){
        term.RE = paste('Point[0-9]+:', 'Group', '.+', sep='')
        term.rows = grep(term.RE, row.names(summary(lme.trk)$tTable))
        data.frame(Point = as.numeric(levels(factor(df$Point))),
                   summary(lme.trk)$tTable[term.rows,])
    } else data.frame()
}
models$tTable = ddply(trk_data, c("Tract", "Hemisphere"), fit_trk_model2)
# Extract a subset like this: subset(models$tTable, Tract=='ILF' & Hemisphere=='L')

# Add weights=varFunc(~I(1/SD^2)) to do weighted fitting

# Model number of streamlines
summary(lm(Streamlines ~ Group + Hemisphere + Tract, data=trk_props_long))

################################################################################
# Permutation testing
# Function to permute groups and return max t-stat distribution
doPerms <- function(df, demog, nPerms=100){
    cat('Permuting data:\n')
    pb   = txtProgressBar(1,nPerms,1, style=3)
    maxT = NULL
    df   = df[,!colnames(df) %in% colnames(demog)[-1]]
    for(i in 1:nPerms){
        # Permute group membership, but respect ID groupings
        perm.labs  = sample(nrow(demog))
        perm.demog = data.frame(ID=demog$ID[perm.labs], demog[,2:ncol(demog)])
        perm.df    = merge(df, perm.demog)
        # Fit models and obtain max t-stat
        tTable        = ddply(perm.df, c("Tract", "Hemisphere"), fit_trk_model2)
        maxT          = c(maxT, max(abs(tTable$t)))
        setTxtProgressBar(pb,i)
    }
    close(pb)
    maxT
}

# Only correct for multiple comparisons across the tracts with a significant
# Point:Group interaction F-test. Set the corrected p-values of other tracts to 1
Ftests          = subset(models$anova, Term=='Point:Group')
sigFtests       = which(Ftests$p.value<0.05)
sigF.Tract      = Ftests$Tract[sigFtests]
sigF.Hemisphere = Ftests$Hemisphere[sigFtests]
sigF.trk_data   = merge(trk_data, data.frame(Tract=sigF.Tract, Hemisphere=sigF.Hemisphere))

maxT = doPerms(sigF.trk_data, demog, nPerms)

# Execute in a terminal to utilize multiple cores and do lots of perms!
#library(multicore)
#doPerm <- function(iPerm, df, demog) {
#   df         = df[,!colnames(df) %in% colnames(demog)[-1]]
#   perm.labs  = sample(nrow(demog))
#   perm.demog = data.frame(ID=demog$ID[perm.labs], demog[,2:ncol(demog)])
#   perm.df    = merge(df, perm.demog)
#   tTable     = ddply(perm.df, c("Tract", "Hemisphere"), fit_trk_model2)
#   max(abs(tTable$t))
#}
#maxT = as.numeric(mclapply(1:10000, doPerm, df=sigF.trk_data, demog=demog, mc.cores=7))
#save(maxT, file=file.path(exptDir, 'maxT.Rdata'))

#load(file.path(exptDir, 'maxT.Rdata'))

# Function to adjust point:group p-values according to where a given t-stat lies
# along the max t-stat distribution
getP <- function(testT, maxT){
    # If you do 1000 permutations and the result from the real dataset is the
    # most extreme, then the empirical p-value is 1/1001
    (1+length(maxT[abs(testT) < maxT]))/(length(maxT)+1)
}
maxT.no.na = maxT[is.na(maxT)==F]
crit       = floor((1-thresh)*length(maxT.no.na))
models$tTable$p.val.adj = sapply(models$tTable$t, getP, maxT.no.na)
models$tTable$p.val.adj[-as.numeric(row.names(subset(models$tTable, Tract==sigF.Tract & Hemisphere==sigF.Hemisphere)))] = 1

# Add a Position column for easier plotting
models$tTable = ddply(models$tTable, c("Tract", "Hemisphere"),
                      transform, Position = (as.numeric(Point)-1) * 100/(max(as.numeric(Point))-1))

################################################################################
# Figures
########
# Plot the distribution of the max t-stat across all points/tracts under
# the null hypothesis of no group effect
dev.new(width=4, height=3)
p1 = ggplot(data.frame(maxT.no.na), aes(x=maxT.no.na, y=..density..))
p1 + geom_histogram(binwidth=0.1) + geom_density(color='blue') + geom_vline(xint = sort(maxT.no.na)[crit], color='red') + xlim(0,5) + ylim(0,1.75) + labs(x=expression(paste('Empirical max ', group('|', italic(t), '|'))), y='Density')

########
# Plot # of streamlines by Group, tract, and hemisphere
p2 = ggplot(trk_props_long, aes(x=Group, y=Streamlines, fill=Group)) + geom_bar(stat='summary', fun.y=mean) + facet_grid(Tract~Hemisphere)

# Set colors
brew_fill = scale_fill_manual(values=rev(brewer.pal(2, 'Set1')[1:2]))

# Generate error bars (pointwise 95% CIs)
stderr   <- function(x) sqrt(var(x)/length(x))
get_bars <- function(y.in){
    data.frame(y=mean(y.in), ymax=mean(y.in)+1.96*stderr(y.in), ymin=mean(y.in)-1.96*stderr(y.in))
}
error_bars = geom_errorbar(stat='summary', fun.data=get_bars, size=0.25, width=0.25)

# Set scales and coordinate systems
set_coords = c(scale_y_continuous(breaks=100*c(0:3)), scale_x_discrete('', breaks=''))

# Add points for individual subjects
jitter  = position_jitter(width=0.2)
sub_pts = geom_point(fill=NA, position=jitter, alpha=0.3)

# Options
set_opts = ylab("Streamlines")

# Make final plot
dev.new(width=3.5, height=4.9)
p2 + set_coords + brew_fill + sub_pts + error_bars + set_opts

########
# Plot FA vs. position, conditioned on hemisphere, tract, and group
p3        = qplot(Position, FA, group=ID, colour=Group, size=Streamlines, alpha=I(0.3), data=trk_data, facets = Tract~Hemisphere, geom='line', xlab='Position along tract (%)') + scale_size(to=c(0.25,3))

# Set colors
brew_cols = scale_colour_manual(values=rev(brewer.pal(2, 'Set1')[1:2]))

# Add group means
means_smooth = stat_smooth(aes(ymax=..y..+1.96*..se.., ymin=..y..-1.96*..se.., group=Group), alpha=0.8, span=0.5)
means        = stat_summary(fun.y=mean, geom='line', size=0.6, aes(group=Group))

# Add an asterisk if there is a significant main group effect
anova_grp  = subset(models$anova, Term=="Group")
Caption    = rep('', nrow(anova_grp))
Caption[anova_grp$p.value<thresh] = '*'
anova_grp  = data.frame(anova_grp, Caption)
grp_effect = geom_text(aes(x=0, y=0.2, label=Caption, group=NULL, size=NULL), data=anova_grp, colour='black')

# If the F-test across the Point:Group terms in a panel is significant, plot a bar at the bottom to indicate which pointwise t-tests are significant
get_breaks <- function(df, thresh=0.05){
    sig  = df$p.value < thresh
    dsig = c(diff(sig), 0)
    if(subset(models$anova, Term=="Point:Group" & Tract==df$Tract[1] & Hemisphere==df$Hemisphere[1])$p.value < thresh){
        onpts  = df$Point[dsig==1]  + 0.5
        offpts = df$Point[dsig==-1] + 0.5
        
        # Check for any unclosed segments
        if(dsig[which(dsig != 0)[1]] == -1) {
            onpts = c(1, onpts)
        }
        if(dsig[rev(which(dsig != 0))[1]] == 1) {
            offpts = c(offpts, length(dsig))
        }
        
        data.frame(on = (onpts-1)  * 100/(length(df$Point)-1),
                  off = (offpts-1) * 100/(length(df$Point)-1))
    } else data.frame(on=0, off=0)
}
break_list = ddply(models$tTable, c("Tract", "Hemisphere"), get_breaks, thresh=thresh)
sig_bars   = geom_segment(aes(x=on, y=0.2, xend=off, yend=0.2, group=NULL, size=NULL), data=break_list, colour='black')

# Make final plot
dev.new(width=7, height=5)
p3 + brew_cols + means_smooth + grp_effect + sig_bars

########
# Plot p-values
# Generate a version of stats results with finer spacing for highlighting
# significant areas in green in p/t-value plots
lin_interp = function(x, spacing=0.01) {
    approx(1:length(x), x, xout=seq(1,length(x), spacing))$y
}
models$tTable.interp = ddply(models$tTable, c("Tract", "Hemisphere"), colwise(lin_interp, c('Point','t.value', 'p.value', 'Position')))

p4 = ggplot(data=models$tTable.interp, aes(x=Position, y=p.value, color=p.value < 0.05, group=1)) + geom_line() + facet_grid(facets=Tract~Hemisphere) + xlab('Position along tract (%)') + scale_color_manual('p.value < 0.05', values=c('red', 'green'))

# Gray out non-significant areas
grayout_p = annotate('rect', xmin=0, xmax=100, ymin=0.05, ymax=1, alpha=0.25)

# Make final plot
dev.new(width=7, height=5)
p4 + grayout_p + scale_y_log10(limits=c(0.001, 1), breaks=c(0.001, 0.01, 0.1, 1))

########
# Plot t-statistics
# Gray out non-significant areas
tcrits = ddply(models$tTable, c("Tract", "Hemisphere"), summarize, t.crit=qt(0.05/2, DF[1]))
grayout_t = geom_rect(aes(x=1, y=1, ymin=t.crit, ymax=-t.crit), color=NA, xmin=0, xmax=100, alpha=0.25, data=tcrits)

# Make final plot
dev.new(width=7, height=5)
p4 + grayout_t + aes(y=t.value, color=t.value < tcrits$t.crit | t.value > -tcrits$t.crit) + ylim(c(-4, 4)) + geom_hline(yint=0, linetype=3)

################################################################################
# Output statistical results for import into MATLAB and overlay onto mean tract # geometry
write.table(models$tTable, file=file.path(exptDir, 'effects_table.txt'), quote=F, row.names=F)