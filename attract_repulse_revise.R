
library(plyr)
library(data.table)
library(ggplot2)
library(scales)
attract = func

sample <- "smc2_asyn_merge"

saddle_bin <- 40 #number of pixels for each comp

saddle_all <- matrix(0, nrow= saddle_bin*4, ncol= saddle_bin*4) #4 comps together.

#create a table to store chromosome compartment_scores
c_score_all <- matrix(0, nrow=20, ncol=11)

colnames(c_score_all) <- c("chr","a1_vs_a1","a1_vs_a2","a1_vs_b1","a1_vs_c","a2_vs_a2","a2_vs_b1","a2_vs_c","b1_vs_b1","b1_vs_c","c_vs_c")

#read in tracks of histone modification as well as PC1 values for each time point
track <- read.table(file="/Users/haoyuezhang/Desktop/szbl/科研项目/smc2-aid/data_analysis/25k_call_compartments/comp_bigwigaveoverbed/his_modi_all_pc1_all_25kb_bins_with_comp_with_smc3.csv",sep=",",header=T)

track <- track[-nrow(track),] #cut the last row from chrX



for(i in 1:20){
	
	print(i)
#matrix to store saddle plot for each chromosome

	if(i==10){next} #skip chr12
#	if(i==2){next} #skip chr2
	if(i==12){next} #skip chr12
	
	saddle <- matrix(0, nrow= saddle_bin*4, ncol= saddle_bin*4)
	
	input_matrix <- paste("/Volumes/super_adam3/smc2_project/matrix_file/25kb/",sample,"_chr", i, ".txt", sep="")

# read in 25kb o/e matrix for each chromosome and rea
	ma <- read.table(input_matrix, sep="\t", header=F)
	
# read in pc1 and histone modification value for each chromosome
	pc1_4h <- track[track$chr == paste("chr",i,sep=""), 12]
	pc1_asy <- track[track$chr == paste("chr",i,sep=""), 39]
	k9me3 <- track[track$chr == paste("chr",i,sep=""), 11]
	k27me3 <- track[track$chr == paste("chr",i,sep=""), 13]
	k36me3 <- track[track$chr == paste("chr",i,sep=""), 6]
	k27ac <- track[track$chr == paste("chr",i,sep=""), 5]
	comp <- track[track$chr == paste("chr",i,sep=""), 65]
	
	length <- nrow(ma)
	ma <- ma[, -(length+1)]
	shuffle <- sample(1:length) #determine a shuffle to order both matrix and comp
	
	order_x <- c(1:length)
	order_y <- c(1:length)
	
	comp_x <- comp[order_x] #order comp based on choosen x
	comp_y <- comp[order_y]
	pc1_asy_x <- pc1_asy[order_x] 
	pc1_asy_y <- pc1_asy[order_y] 
	ma <- ma[order_x, order_y]
		
	y <- comp_y
	x <- comp_x
	
	
# add track value to each row of the oe matrix for y axis sorting. Note the first y_sort is columne name. 
	ma$y_sort <- y[c(1:length)]
	ma$pc1 <- pc1_asy_y[c(1:length)]

# delete rows that has pc1 =0, which means unsequenced region, also, only focus on A_1, A_2, B_1 and C compartments
	ma <- ma[ma$pc1 !=0,]

# sort each row based on comp in ascending order

	ma <- ma[order(ma$y_sort, decreasing=F),]

# remove the pc1 column from o/e matrix. but first, temperorily store the y_sort columne

	tmp <- ma$y_sort #this is to indicate which bins belong to which compartment in the future collapsing step.
	ma <- ma[, -c((length+1):(length+2))]
	
# add a row of pc1 to label each column. This is for x axis sorting
	ma <- rbind(ma, x[c(1:length)])
	ma <- rbind(ma,pc1_asy_x[c(1:length)])

# remove columns with x =0, which means unsequenced region. also remove b_2 AND d compartments
	ma <- ma[,-which(pc1_asy_x==0)]
	
# sort each column based on x value.

	ma <- ma[,order(ma[(nrow(ma)-1),], decreasing =F)] #note, has to order based on the second to last row. because pc1 is the last row

# remove the the row with x values
	ma <- ma[-nrow(ma),]
	ma <- ma[-nrow(ma),] #do this twice to remove the last row and second to last row

	ma$comp <- tmp # add the tmp column back to ma matrix to indicate comp bins

# collapse every comp into a matrix with saddle_bin of length. This is to average every chromosome together.
# the residual pixels that cannot be equally divided will put at the beginning.

	a_1_l <- nrow(ma[ma$comp=="A_1",])
	a_2_l <- nrow(ma[ma$comp=="A_2",])
	b_1_l <- nrow(ma[ma$comp=="B_1",])
	b_2_l <- nrow(ma[ma$comp=="B_2",]) #note, there's a B2 between b1 and that should be considered !
	c_l <- nrow(ma[ma$comp=="B_3",]) #note, here, C is defined as B_3 in the table

bin_perU_a_1 <- floor(a_1_l/(saddle_bin)) #note that it has to be divided by "saddle_bin-1"
bin_perU_a_2 <- floor(a_2_l/(saddle_bin)) #note that it has to be divided by "saddle_bin-1"
bin_perU_b_1 <- floor(b_1_l/(saddle_bin)) #note that it has to be divided by "saddle_bin-1"
bin_perU_c <- floor(c_l/(saddle_bin)) #note that it has to be divided by "saddle_bin-1"

a_1_bin1_s <- 1
a_1_bin1_e <- bin_perU_a_1
a_2_bin1_s <- a_1_l +1
a_2_bin1_e <- a_1_l + bin_perU_a_2
b_1_bin1_s <- a_1_l + a_2_l +1
b_1_bin1_e <- a_1_l + a_2_l + bin_perU_b_1
c_bin1_s <- a_1_l + a_2_l + b_1_l + b_2_l + 1 #note, had to add b_2_l
c_bin1_e <- a_1_l + a_2_l + b_1_l + b_2_l + bin_perU_c #note, had to add b_2_l



#fill in the averaged value

	for(j in 1: saddle_bin){ #4 sections adding together. therefore, the matrix should be saddle_bin*4. I only need to fill in the upper triangle.
		
		print(j)
		
		for(k in 1: saddle_bin){	 #only filling the upper triangle	
			
			if(j==1 & k==1){
						
			saddle[j,k] <- mean(colMeans(ma[c(a_1_bin1_s: a_1_bin1_e),c(a_1_bin1_s: a_1_bin1_e)], na.rm=T))
			saddle[j,(k+ saddle_bin*1)] <- mean(colMeans(ma[c(a_1_bin1_s: a_1_bin1_e), c(a_2_bin1_s: a_2_bin1_e)], na.rm=T))
			saddle[j,(k+ saddle_bin*2)] <- mean(colMeans(ma[c(a_1_bin1_s: a_1_bin1_e), c(b_1_bin1_s: b_1_bin1_e)], na.rm=T))
			saddle[j,(k+ saddle_bin*3)] <- mean(colMeans(ma[c(a_1_bin1_s: a_1_bin1_e), c(c_bin1_s: c_bin1_e)], na.rm=T))
			
			saddle[(j+saddle_bin*1),(k+saddle_bin*1)] <- mean(colMeans(ma[c(a_2_bin1_s: a_2_bin1_e),c(a_2_bin1_s: a_2_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*1),(k+saddle_bin*2)] <- mean(colMeans(ma[c(a_2_bin1_s: a_2_bin1_e), c(b_1_bin1_s: b_1_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*1),(k+saddle_bin*3)] <- mean(colMeans(ma[c(a_2_bin1_s: a_2_bin1_e),c(c_bin1_s: c_bin1_e)], na.rm=T))

			saddle[(j+saddle_bin*2),(k+saddle_bin*2)] <- mean(colMeans(ma[c(b_1_bin1_s: b_1_bin1_e), c(b_1_bin1_s: b_1_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*2),(k+saddle_bin*3)] <- mean(colMeans(ma[c(b_1_bin1_s: b_1_bin1_e),c(c_bin1_s: c_bin1_e)], na.rm=T))
			
			saddle[(j+saddle_bin*3),(k+saddle_bin*3)] <- mean(colMeans(ma[c(c_bin1_s: c_bin1_e),c(c_bin1_s: c_bin1_e)], na.rm=T))
			
			saddle[(j+saddle_bin*1),k] <- mean(colMeans(ma[c(a_2_bin1_s: a_2_bin1_e),c(a_1_bin1_s: a_1_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*2),k] <- mean(colMeans(ma[c(b_1_bin1_s: b_1_bin1_e),c(a_1_bin1_s: a_1_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*3),k] <- mean(colMeans(ma[c(c_bin1_s: c_bin1_e),c(a_1_bin1_s: a_1_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*2),(k+saddle_bin*1)] <- mean(colMeans(ma[c(b_1_bin1_s: b_1_bin1_e),c(a_2_bin1_s: a_2_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*3),(k+saddle_bin*1)] <- mean(colMeans(ma[c(c_bin1_s: c_bin1_e),c(a_2_bin1_s: a_2_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*3),(k+saddle_bin*2)] <- mean(colMeans(ma[c(c_bin1_s: c_bin1_e),c(b_1_bin1_s: b_1_bin1_e)], na.rm=T))
					
			
			}
			
			if(j==1 & k !=1){
			a1_s <- (k-2)*bin_perU_a_1 + a_1_bin1_e +1
			a1_e <- (k-1)* bin_perU_a_1 + a_1_bin1_e
			a2_s <- (k-2)*bin_perU_a_2 + a_2_bin1_e +1
			a2_e <- (k-1)* bin_perU_a_2 + a_2_bin1_e
			b1_s <- (k-2)*bin_perU_b_1 + b_1_bin1_e +1
			b1_e <- (k-1)* bin_perU_b_1 + b_1_bin1_e
			c_s <- (k-2)*bin_perU_c + c_bin1_e +1
			c_e <- (k-1)* bin_perU_c + c_bin1_e

			saddle[j,k] <- mean(colMeans(ma[c(a_1_bin1_s: a_1_bin1_e),c(a1_s:a1_e)], na.rm=T))
			saddle[j,(k+saddle_bin*1)] <- mean(colMeans(ma[c(a_1_bin1_s: a_1_bin1_e), c(a2_s: a2_e)], na.rm=T))
			saddle[j,(k+saddle_bin*2)] <- mean(colMeans(ma[c(a_1_bin1_s: a_1_bin1_e), c(b1_s: b1_e)], na.rm=T))
			saddle[j,(k+saddle_bin*3)] <- mean(colMeans(ma[c(a_1_bin1_s: a_1_bin1_e), c(c_s: c_e)], na.rm=T))

			saddle[(j+saddle_bin*1),(k+saddle_bin*1)] <- mean(colMeans(ma[c(a_2_bin1_s: a_2_bin1_e),c(a2_s: a2_e)], na.rm=T))
			saddle[(j+saddle_bin*1),(k+saddle_bin*2)] <- mean(colMeans(ma[c(a_2_bin1_s: a_2_bin1_e), c(b1_s: b1_e)], na.rm=T))
			saddle[(j+saddle_bin*1),(k+saddle_bin*3)] <- mean(colMeans(ma[c(a_2_bin1_s: a_2_bin1_e),c(c_s: c_e)], na.rm=T))

			saddle[(j+saddle_bin*2),(k+saddle_bin*2)] <- mean(colMeans(ma[c(b_1_bin1_s: b_1_bin1_e), c(b1_s: b1_e)], na.rm=T))
			saddle[(j+saddle_bin*2),(k+saddle_bin*3)] <- mean(colMeans(ma[c(b_1_bin1_s: b_1_bin1_e),c(c_s: c_e)], na.rm=T))
			
			saddle[(j+saddle_bin*3),(k+saddle_bin*3)] <- mean(colMeans(ma[c(c_bin1_s: c_bin1_e),c(c_s: c_e)], na.rm=T))
			
			saddle[(j+saddle_bin*1),k] <- mean(colMeans(ma[c(a_2_bin1_s: a_2_bin1_e),c(a1_s:a1_e)], na.rm=T))
			saddle[(j+saddle_bin*2),k] <- mean(colMeans(ma[c(b_1_bin1_s: b_1_bin1_e),c(a1_s:a1_e)], na.rm=T))
			saddle[(j+saddle_bin*3),k] <- mean(colMeans(ma[c(c_bin1_s: c_bin1_e),c(a1_s:a1_e)], na.rm=T))
			saddle[(j+saddle_bin*2),(k+saddle_bin*1)] <- mean(colMeans(ma[c(b_1_bin1_s: b_1_bin1_e),c(a2_s: a2_e)], na.rm=T))
			saddle[(j+saddle_bin*3),(k+saddle_bin*1)] <- mean(colMeans(ma[c(c_bin1_s: c_bin1_e),c(a2_s: a2_e)], na.rm=T))
			saddle[(j+saddle_bin*3),(k+saddle_bin*2)] <- mean(colMeans(ma[c(c_bin1_s: c_bin1_e),c(b1_s: b1_e)], na.rm=T))
							
			}



			if(j!=1 & k ==1){
			a1_s_h <- (j-2)*bin_perU_a_1 + a_1_bin1_e +1
			a1_e_h <- (j-1)* bin_perU_a_1 + a_1_bin1_e
			a2_s_h <- (j-2)*bin_perU_a_2 + a_2_bin1_e +1
			a2_e_h <- (j-1)* bin_perU_a_2 + a_2_bin1_e
			b1_s_h <- (j-2)*bin_perU_b_1 + b_1_bin1_e +1
			b1_e_h <- (j-1)* bin_perU_b_1 + b_1_bin1_e
			c_s_h <- (j-2)*bin_perU_c + c_bin1_e +1
			c_e_h <- (j-1)* bin_perU_c + c_bin1_e

			saddle[j,k] <- mean(colMeans(ma[c(a1_s_h: a1_e_h),c(a_1_bin1_s: a_1_bin1_e)], na.rm=T))
			saddle[j,(k+saddle_bin*1)] <- mean(colMeans(ma[c(a1_s_h: a1_e_h), c(a_2_bin1_s: a_2_bin1_e)], na.rm=T))
			saddle[j,(k+saddle_bin*2)] <- mean(colMeans(ma[c(a1_s_h: a1_e_h), c(b_1_bin1_s: b_1_bin1_e)], na.rm=T))
			saddle[j,(k+saddle_bin*3)] <- mean(colMeans(ma[c(a1_s_h: a1_e_h), c(c_bin1_s: c_bin1_e)], na.rm=T))

			saddle[(j+saddle_bin*1),(k+saddle_bin*1)] <- mean(colMeans(ma[c(a2_s_h: a2_e_h),c(a_2_bin1_s: a_2_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*1),(k+saddle_bin*2)] <- mean(colMeans(ma[c(a2_s_h: a2_e_h),c(b_1_bin1_s: b_1_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*1),(k+saddle_bin*3)] <- mean(colMeans(ma[c(a2_s_h: a2_e_h),c(c_bin1_s: c_bin1_e)], na.rm=T))

			saddle[(j+saddle_bin*2),(k+saddle_bin*2)] <- mean(colMeans(ma[c(b1_s_h: b1_e_h), c(b_1_bin1_s: b_1_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*2),(k+saddle_bin*3)] <- mean(colMeans(ma[c(b1_s_h: b1_e_h),c(c_bin1_s: c_bin1_e)], na.rm=T))
			
			saddle[(j+saddle_bin*3),(k+saddle_bin*3)] <- mean(colMeans(ma[c(c_s_h: c_e_h),c(c_bin1_s: c_bin1_e)], na.rm=T))
			
			saddle[(j+saddle_bin*1),k] <- mean(colMeans(ma[c(a2_s_h: a2_e_h),c(a_1_bin1_s: a_1_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*2),k] <- mean(colMeans(ma[c(b1_s_h: b1_e_h),c(a_1_bin1_s: a_1_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*3),k] <- mean(colMeans(ma[c(c_s_h: c_e_h),c(a_1_bin1_s: a_1_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*2),(k+saddle_bin*1)] <- mean(colMeans(ma[c(b1_s_h: b1_e_h),c(a_1_bin1_s: a_1_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*3),(k+saddle_bin*1)] <- mean(colMeans(ma[c(c_s_h: c_e_h),c(a_2_bin1_s: a_2_bin1_e)], na.rm=T))
			saddle[(j+saddle_bin*3),(k+saddle_bin*2)] <- mean(colMeans(ma[c(c_s_h: c_e_h),c(b_1_bin1_s: b_1_bin1_e)], na.rm=T))
							
			}


			if(j !=1 & k !=1){
			
			a1_s <- (k-2)*bin_perU_a_1 + a_1_bin1_e +1
			a1_e <- (k-1)* bin_perU_a_1 + a_1_bin1_e
			a2_s <- (k-2)*bin_perU_a_2 + a_2_bin1_e +1
			a2_e <- (k-1)* bin_perU_a_2 + a_2_bin1_e
			b1_s <- (k-2)*bin_perU_b_1 + b_1_bin1_e +1
			b1_e <- (k-1)* bin_perU_b_1 + b_1_bin1_e
			c_s <- (k-2)*bin_perU_c + c_bin1_e +1
			c_e <- (k-1)* bin_perU_c + c_bin1_e

			a1_s_h <- (j-2)*bin_perU_a_1 + a_1_bin1_e +1
			a1_e_h <- (j-1)* bin_perU_a_1 + a_1_bin1_e
			a2_s_h <- (j-2)*bin_perU_a_2 + a_2_bin1_e +1
			a2_e_h <- (j-1)* bin_perU_a_2 + a_2_bin1_e
			b1_s_h <- (j-2)*bin_perU_b_1 + b_1_bin1_e +1
			b1_e_h <- (j-1)* bin_perU_b_1 + b_1_bin1_e
			c_s_h <- (j-2)*bin_perU_c + c_bin1_e +1
			c_e_h <- (j-1)* bin_perU_c + c_bin1_e

			saddle[j,k] <- mean(colMeans(ma[c(a1_s_h: a1_e_h),c(a1_s:a1_e)], na.rm=T))
			saddle[j,(k+saddle_bin*1)] <- mean(colMeans(ma[c(a1_s_h: a1_e_h), c(a2_s: a2_e)], na.rm=T))
			saddle[j,(k+saddle_bin*2)] <- mean(colMeans(ma[c(a1_s_h: a1_e_h), c(b1_s: b1_e)], na.rm=T))
			saddle[j,(k+saddle_bin*3)] <- mean(colMeans(ma[c(a1_s_h: a1_e_h), c(c_s: c_e)], na.rm=T))
			
			saddle[(j+saddle_bin*1),(k+saddle_bin*1)] <- mean(colMeans(ma[c(a2_s_h: a2_e_h),c(a2_s: a2_e)], na.rm=T))
			saddle[(j+saddle_bin*1),(k+saddle_bin*2)] <- mean(colMeans(ma[c(a2_s_h: a2_e_h), c(b1_s: b1_e)], na.rm=T))
			saddle[(j+saddle_bin*1),(k+saddle_bin*3)] <- mean(colMeans(ma[c(a2_s_h: a2_e_h),c(c_s: c_e)], na.rm=T))

			saddle[(j+saddle_bin*2),(k+saddle_bin*2)] <- mean(colMeans(ma[c(b1_s_h: b1_e_h), c(b1_s: b1_e)], na.rm=T))
			saddle[(j+saddle_bin*2),(k+saddle_bin*3)] <- mean(colMeans(ma[c(b1_s_h: b1_e_h),c(c_s: c_e)], na.rm=T))
			
			saddle[(j+saddle_bin*3),(k+saddle_bin*3)] <- mean(colMeans(ma[c(c_s_h: c_e_h),c(c_s: c_e)], na.rm=T))
			
			saddle[(j+saddle_bin*1),k] <- mean(colMeans(ma[c(a2_s_h: a2_e_h),c(a1_s:a1_e)], na.rm=T))
			saddle[(j+saddle_bin*2),k] <- mean(colMeans(ma[c(b1_s_h: b1_e_h),c(a1_s:a1_e)], na.rm=T))
			saddle[(j+saddle_bin*3),k] <- mean(colMeans(ma[c(c_s_h: c_e_h),c(a1_s:a1_e)], na.rm=T))
			saddle[(j+saddle_bin*2),(k+saddle_bin*1)] <- mean(colMeans(ma[c(b1_s_h: b1_e_h),c(a2_s: a2_e)], na.rm=T))
			saddle[(j+saddle_bin*3),(k+saddle_bin*1)] <- mean(colMeans(ma[c(c_s_h: c_e_h),c(a2_s: a2_e)], na.rm=T))
			saddle[(j+saddle_bin*3),(k+saddle_bin*2)] <- mean(colMeans(ma[c(c_s_h: c_e_h),c(b1_s: b1_e)], na.rm=T))



			}			
			
		}
		
}

# compute the interactions of each comp among chr

	c_score_all[i,1] <- paste("chr",i,sep="")
	c_score_all[i,2] <- mean(saddle[c(1: saddle_bin),c(1: saddle_bin)])
	c_score_all[i,3] <- mean(saddle[c(1: saddle_bin),c((saddle_bin+1): (saddle_bin*2))])
	c_score_all[i,4] <- mean(saddle[c(1: saddle_bin),c((saddle_bin*2+1):(saddle_bin*3))])
	c_score_all[i,5] <- mean(saddle[c(1: saddle_bin),c((saddle_bin*3+1): (saddle_bin*4))])
	c_score_all[i,6] <- mean(saddle[c((saddle_bin+1): (saddle_bin*2)),c((saddle_bin+1): (saddle_bin*2))])
	c_score_all[i,7] <- mean(saddle[c((saddle_bin+1): (saddle_bin*2)),c((saddle_bin*2+1):(saddle_bin*3))])
	c_score_all[i,8] <- mean(saddle[c((saddle_bin+1): (saddle_bin*2)),c((saddle_bin*3+1):(saddle_bin*4))])
	c_score_all[i,9] <- mean(saddle[c((saddle_bin*2+1): (saddle_bin*3)),c((saddle_bin*2+1):(saddle_bin*3))])
	c_score_all[i,10] <- mean(saddle[c((saddle_bin*2+1): (saddle_bin*3)),c((saddle_bin*3+1):(saddle_bin*4))])
	c_score_all[i,11] <- mean(saddle[c((saddle_bin*3+1): (saddle_bin*4)),c((saddle_bin*3+1):(saddle_bin*4))])

	
	
#	c_score_all[i,3] <- lr_corner
#	c_score_all[i,4] <- ur_corner
#	c_score_all[i,5] <- ll_corner
#	c_score_all[i,6] <- mean(ul_corner,lr_corner)/mean(ur_corner,ll_corner)
	
	saddle_all <- saddle_all + saddle

}

saddle_all <- saddle_all/18 

write.table(saddle_all, file=paste("/Users/haoyuezhang/Desktop/szbl/科研项目/smc2-aid/data_analysis/25k_call_compartments/comp_plots/interpolation/",sample,"_comp_original_order_18chr",".txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F)

write.table(c_score_all, file=paste("/Users/haoyuezhang/Desktop/szbl/科研项目/smc2-aid/data_analysis/25k_call_compartments/comp_plots/interpolation/",sample,"_comp_interactions_original_order_18chr",".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

#generate genomewide log2 transformed saddle plot for generation

saddle_all_plot <- log2(saddle_all)

colnames(saddle_all_plot) <- NULL
rownames(saddle_all_plot) <- NULL

# draw plot

saddle_all_plot.melt <- melt(saddle_all_plot)

## plot matrix
png(file=paste("/Users/haoyuezhang/Desktop/szbl/科研项目/smc2-aid/data_analysis/25k_call_compartments/comp_plots/interpolation/",sample,"_original_order_comp_18chr",".png",sep=""))
ggplot(saddle_all_plot.melt) + geom_tile(aes(Var1,-Var2,fill=value)) + scale_fill_gradient2(low="#006CB9", high="#B40003", mid="white", midpoint=0, limits=c(-2,2), oob=squish) + theme_minimal() + theme(axis.text=element_blank()) + coord_fixed() + xlab("") + ylab("") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
dev.off()

