

sub convert_IUPAC_string_to_nucleotides()
	{
		my ($user_repeat,$ref_lib_file)=@_;
		
		my $return_repeat_string="";
		
		#if($user_repeat ne "Not_Specified")
		#	{
				$user_repeat=uc($user_repeat);
				
				my @arr_user_repeat=split('',$user_repeat);
				my @user_repeat_parts;
					
				foreach my $base(@arr_user_repeat)
					{
						my @tmp_arr;
						
						if($base eq "N")
							{
								push(@tmp_arr,'A');
								push(@tmp_arr,'C');
								push(@tmp_arr,'G');
								push(@tmp_arr,'T');
							}
						elsif($base eq "R")
							{
								push(@tmp_arr,'A');
								push(@tmp_arr,'G');
							}
						elsif($base eq "Y")
							{
								push(@tmp_arr,'C');
								push(@tmp_arr,'T');
							}
						elsif($base eq "S")
							{
								push(@tmp_arr,'G');
								push(@tmp_arr,'C');
							}
						elsif($base eq "W")
							{
								push(@tmp_arr,'A');
								push(@tmp_arr,'T');
							}
						elsif($base eq "K")
							{
								push(@tmp_arr,'G');
								push(@tmp_arr,'T');
							}
						elsif($base eq "M")
							{
								push(@tmp_arr,'A');
								push(@tmp_arr,'C');
							}
						elsif($base eq "B")
							{
								push(@tmp_arr,'C');
								push(@tmp_arr,'G');
								push(@tmp_arr,'T');
							}
						elsif($base eq "D")
							{
								push(@tmp_arr,'A');
								push(@tmp_arr,'G');
								push(@tmp_arr,'T');
							}
						elsif($base eq "H")
							{
								push(@tmp_arr,'A');
								push(@tmp_arr,'C');
								push(@tmp_arr,'T');
							}
						elsif($base eq "V")
							{
								push(@tmp_arr,'A');
								push(@tmp_arr,'C');
								push(@tmp_arr,'G');
							}
						if($base eq ".")
							{
								push(@tmp_arr,'-');		
							}
						if($base eq "-")
							{
								push(@tmp_arr,'-');		
							}	
						else{
								push(@tmp_arr,$base);
							}		
						
						
						#------------------------
						#my $key_index=0;
						
						if($#user_repeat_parts<0)
							{
								foreach my $nt(@tmp_arr)
									{
										push(@user_repeat_parts,$nt);
									}
							}
						else{
								my @tmp_arr2;
								for my $i(0..$#user_repeat_parts)
									{
										my $old_pam_part=shift(@user_repeat_parts);
										#my $old_pam_part=$user_repeat_parts[$i];
										
										foreach my $nt(@tmp_arr)
											{
												my $new_pam_part=$old_pam_part.$nt;												
												push(@tmp_arr2,$new_pam_part);
												
												#print "$new_pam_part<br>";
											}		
									}
								@user_repeat_parts=@tmp_arr2;	
							}

							
					}
				
				
				
				#my $new_user_repeats="";
				foreach my $pam(@user_repeat_parts)
					{
						#$new_user_repeats=$new_user_repeats.$pam.",";
						$return_repeat_string=$return_repeat_string.$pam.",";
						#print "$pam<br>";	
						
									
					}
				#system("echo '$return_repeat_string\tRETURN_REPEAT' >>tmp/$ref_lib_file");	
				
				#print "\$user_repeat=$user_repeat<br>\$new_user_repeats=$new_user_repeats<br>";
				
				#$hash_of_pams_3p{'CRISPR-user'}=$new_user_repeats;
		#	}
		
		
		#---------------------------- return return_repeat_string ------------------------------------------
		
		return($return_repeat_string);
	}



sub process_gbk_file()
	{
		my($input_gbk_file,$output_gene_position_file)=@_;
		
		my @arr_tmp_1=split("/",$input_gbk_file);
		
		my $accession=$arr_tmp_1[$#arr_tmp_1]; $accession=~s/\.\S+$//;
		
		#print "Accession : $accession [$input_gbk_file] \$input_gbk_file=$input_gbk_file \n";
		
		#-------------------------------------------------------------------------------
		open(RD,"tmp/$input_gbk_file") or print "$!";
		my @arr_gbk_file=<RD>;
		close(RD);
		
		#------  extract fasta seq ------------------------------------------------------
		#system("seqret -sequence $input_gbk_file -outseq $output_fasta_file");
		#------ extract the gene positions ----------------------------------------------
		open(TAB,">>tmp/$output_gene_position_file") or print "";
		flock(TAB,2);
		
		my @arr_gene_info_line=`grep -nH -w 'CDS' tmp/$input_gbk_file >&1`;
		#print "grep -nH -w '     CDS' tmp/$input_gbk_file\n @arr_gene_info_line";exit;
		
		foreach my $gene_info_line(@arr_gene_info_line)
			{
				if($gene_info_line=~/(\d+\.\.\d+)/)
					{
						my($start,$stop)=split('\.\.',$1);
						
						print TAB "$accession\tCDS\t$start\t$stop\tNA\tNA\n";
					}
			}
		
		
		
		#-------- extract Cas genes ------------------------------------------------------
		my @arr_lines=`grep -nH 'CRISPR' tmp/$input_gbk_file >&1`;
		
		my %hash_of_crispr_regions;
		foreach my $rec_line(@arr_lines)
			{
				#print "\n\n$rec_line";
				my($t_file,$line_number,$tmp_1)=split(':',$rec_line);
				
				my $cas_gene_start_stop;
				my $cas_gene_start;
				my $cas_gene_stop;
				my $cas_gene_found=0;
				my $i=$line_number-1;
				
				my $skip_this_rec=0;
				while($i>0)
					{
						my $current_line=$arr_gbk_file[$i];chomp $current_line; $current_line=~s/\r+$//;

						if($current_line=~/^     \S+\s+/)
							{
								#print "\t$i :$current_line\n";								
								#---------------------------------
								$current_line=~s/>//g;$current_line=~s/<//g;
								if($current_line=~/(\d+\.\.\d+)/)
									{
										$cas_gene_start_stop=$1;
										
										($cas_gene_start,$cas_gene_stop)=split('\.\.',$1);	
																			
										#print "\tCas gene start-stop=$cas_gene_start-$cas_gene_stop\n";
										$cas_gene_found=1;
									}
								else{
										$skip_this_rec=1;  #wrong rec, should not proceed further.
										last;
									}	
							}
							
							
						if($cas_gene_found!=0){last;}		
						$i--;
					}
				if($skip_this_rec==1){next;}
				if(not $cas_gene_start or not $cas_gene_stop){next;}
					
				#--------- now get the Cas family -----------------------------------	
				my $cas_gene_family_found=0;
				my $cas_det_line="";	
							
				my $j=$line_number-1; #---- as grep line number starts with 1 -------
				while($j<$#arr_gbk_file)
					{
						my $current_line=$arr_gbk_file[$j];chomp $current_line; $current_line=~s/^\s+//; $current_line=~s/\r+//g;
						my $next_line=$arr_gbk_file[$j+1];$next_line=~s/^\s+//;
						
						$cas_det_line=$cas_det_line.$current_line;
						
						if($next_line=~/^\//)
							{
								#$cas_det_line=$cas_det_line.$current_line;
								last;
							}
						elsif($next_line=~/^ORIGIN/)
							{
								last;
							}	
	
	
							
							
						if($cas_gene_family_found!=0){last;}		
						$j++;
					}
					
				$cas_det_line=~s/\/note=//g;$cas_det_line=~s/\"//g;$cas_det_line=~s/\'//g;		
				#print "\t$cas_det_line\n";
				
				#---- now extract the family type ----------
				my $cas_family="NA";
				my $cas_family_type="NA";
				my @arr_family_det=split(';',$cas_det_line);
				if($arr_family_det[0]=~/(C\D\D\d+)/)
					{
						$cas_family=$1; $cas_family=~s/\.\./_/g;  $cas_family=~s/^\s+//; $cas_family=~s/^,//;
						
						if($cas_family eq "")
							{
								$cas_family="NA";
							}
						#print "\tFamily:$cas_family\n";
					}	
				if(defined $arr_family_det[1] and $arr_family_det[1]=~/$cas_family(\S+)/)
					{
						$cas_family_type=$1; $cas_family_type=~s/^_//; $cas_family_type=~s/_/,/g;  $cas_family_type=~s/^\s+//; $cas_family_type=~s/^,//;
						
						if($cas_family_type eq "")
							{
								$cas_family_type="NA";
							}
						#print "\tType:$cas_family_type\n";
					}
				#print "\n";
				
				#------- write the CRISPR position, Family and type to the table
				
				#print TAB "$accession\tCRISPR\t$cas_gene_start\t$cas_gene_stop\t$cas_family\t$cas_family_type\n";
				
				#----------- now store the record in %hash_of_crispr_regions ------------------------------------
				if(defined $hash_of_crispr_regions{$cas_gene_start}{$cas_gene_stop})
					{
						my($fam,$type)=split('\t',$hash_of_crispr_regions{$cas_gene_start}{$cas_gene_stop});
						
						my $better_fam="NA"; my $better_type="NA";
						#------- check existing record ---------------------------------------------------------
						
						if($cas_family eq "NA" and $fam ne "NA"){$better_fam=$fam;}
						elsif($fam eq "NA" and $cas_family ne "NA"){$better_fam=$cas_family;}
						elsif($fam eq "NA" and $cas_family eq "NA"){$better_fam="NA";}
						
						if($cas_family_type eq "NA" and $type ne "NA"){$better_type=$type;}
						elsif($type eq "NA" and $cas_family_type ne "NA"){$better_type=$cas_family_type;}
						elsif($type eq "NA" and $cas_family_type eq "NA"){$better_type="NA";}
						
						$hash_of_crispr_regions{$cas_gene_start}{$cas_gene_stop}="$better_fam\t$better_type";
					}
				else{
						$hash_of_crispr_regions{$cas_gene_start}{$cas_gene_stop}="$cas_family\t$cas_family_type";
					}	
			}
			
		#---------------------- now print the CRISPR gene positions --------------------------------------------	
		foreach my $start_p(sort{$a<=>$b} keys %hash_of_crispr_regions)
			{
				foreach my $stop_p(sort{$a<=>$b} keys %{$hash_of_crispr_regions{$start_p}})
					{					
						print TAB "$accession\tCRISPR\t$start_p\t$stop_p\t$hash_of_crispr_regions{$start_p}{$stop_p}\n";
					}	
			}	
			
		close(TAB);
		return 1;
	}



sub convert_crt_arrays_to_pilercr_and_append()
	{
		my($crt_output_file,$pilercr_output_file)=@_;
		
		#----- open the PILER-CR output file and get the array index --------
		open(RD,"tmp/$pilercr_output_file") or print "$!";
		my @arr_pilercr=<RD>;
		close(RD);
		my $highest_array_index=0;
		my %pilercr_array_ranges;
		for(my $i=0;$i<$#arr_pilercr;$i++)
			{
				my $line=$arr_pilercr[$i];
				if($line=~/^Array/ and $arr_pilercr[$i+1]=~/^>/)
					{
						if($line=~/^Array (\d+)/)
							{
								if($1>$highest_array_index){$highest_array_index=$1;}
							}
						#------ now get the start and stop ---------------------------
						my $j=$i+5;
						my $array_start;my $array_stop;
						while($arr_pilercr[$j]!~/====/)
							{
								if($j==$i+5)
									{
										my $first_line=$arr_pilercr[$j];$first_line=~s/^\s+//;
										if($first_line=~/^(\d+) /){$array_start=$1;}
									}
								if($arr_pilercr[$j+1]=~/====/)
									{
										#print "$array_start-$arr_pilercr[$j+1]<br>";	#exit;
										my $last_line=$arr_pilercr[$j];$last_line=~s/^\s+//;$last_line=~s/\s+/ /g;
										
										if($last_line=~/^(\d+) (\d+)/){$array_stop=$1+$2;}
										if($array_stop>$array_start)
										{
											my $range="$array_start-$array_stop";#print "$array_start-$array_stop<br>";
											$pilercr_array_ranges{$range}=1;
										}
										last;
									}
								$j++;		
							}	
						
						#-----------------------------------------------------------------	
						$i=$j;	
					}
			}
			
		#print "\$highest_array_index=$highest_array_index<BR>";#exit;	
		
		#---- now read the arrays predicted by CRT and convert them to piler-CR format --------------
		open(CRT,"tmp/$crt_output_file") or print "$!";
		my @arr_crt=<CRT>;
		close(CRT);
		
		my $accession="";
		for(my $i=0;$i<=$#arr_crt;$i++)
			{
				my $line=$arr_crt[$i]; chomp $line; $line=~s/\r+//g;
				
				if($line=~/ORGANISM:  (\S+)/)
					{
						$accession=$1;
						$highest_array_index=0;
						#print "\$accession=$accession\t\$i=$i\n";
					}
					
				elsif($arr_crt[$i]=~/^CRISPR/ and $arr_crt[$i]=~/Range:/)
					{
						my $range;my $range_center;

						if($arr_crt[$i]=~/Range: (\d+) - (\d+)/)
							{
								$range=$1."-".$2;
								$range_center=$1+int(($2-$1)/2);
							}	
						#--------- check if the array center belong to any piler-cr array already ----
						my $array_exist=0;
						if(defined $range)
							{
								foreach my $pilercr_range(keys %pilercr_array_ranges)
									{
										#print "$pilercr_range\t$range_center<br>";
										my($r_start,$r_stop)=split('-',$pilercr_range);
										if($range_center>=$r_start and $range_center<=$r_stop){$array_exist=1;last;}
									}
							}
						if($array_exist==1){next;}	#exit;	
						#-------- extract all the repeats, and get the model repeat ------------------
						my $j=$i+3;
						my $c_line="";
						my @arr_repeats;
						my %hash_of_repeats;
						while($arr_crt[$j]!~/^---/)
							{
								my @arr_line=split('\t+',$arr_crt[$j]);
								
								push(@arr_repeats,$arr_line[1]);
								
								if($hash_of_repeats{$arr_line[1]})
									{
										$hash_of_repeats{$arr_line[1]}=$hash_of_repeats{$arr_line[1]}+1;
									}
								else{
										$hash_of_repeats{$arr_line[1]}=1;
									}
								
								$j++;
							}
						
						#---------- now get the model repeat ----------------------------------------
						my $model_repeat="";
						my @arr_of_model_repeats;
						my $heighest_repeat_occurrence_score=0;
						foreach my $repeat(sort{$hash_of_repeats{$b}<=>$hash_of_repeats{$a}} keys %hash_of_repeats)
							{
								if(not $arr_of_model_repeats[0])
									{
											#push(@arr_of_model_repeats,$repeat);
										$arr_of_model_repeats[0]=$repeat;
										$heighest_repeat_occurrence_score=$hash_of_repeats{$repeat};
										$model_repeat=$repeat;
									}
								elsif($hash_of_repeats{$repeat}>=$heighest_repeat_occurrence_score)
									{
										push(@arr_of_model_repeats,$repeat);
										$model_repeat=$repeat;
									}
							}
							
						#--------- check if there is just one repeat in the array
						
						if($#arr_of_model_repeats>0)
							{
								#----- now select the best model repeat by scoring the array
								my %tmp_hash1;
								foreach my $repeat(@arr_of_model_repeats)
									{
										my $array_degeneracy_score=&get_array_degeneracy_score_pcr($repeat,\@arr_repeats);
										
										$tmp_hash1{$repeat}=$array_degeneracy_score;
									}
								
								foreach my $repeat(sort{$tmp_hash1{$b}<=>$tmp_hash1{$a}}keys %tmp_hash1)
									{
										#print "$repeat\t$tmp_hash1{$repeat}\n";
										
										$model_repeat=$repeat;
										last;
									}												
							}
						else{
								$model_repeat=$arr_of_model_repeats[0];
							}		
						#print "$model_repeat";exit;	
						
						#----- now that the model repeat is found, convert the whole array on the go
						open(WR,">>tmp/$pilercr_output_file") or print "$!";
						$highest_array_index++;
						print WR "\n";
						print WR "Array $highest_array_index\n";
						print WR ">$accession\n";
						print WR "\n";
						print WR "       Pos  Repeat     %id  Spacer  Left flank    Repeat                               Spacer\n";
						
						my $tmp_equal_sign="";
						for(my $m=0;$m<length($model_repeat);$m++)
							{
								$tmp_equal_sign=$tmp_equal_sign."=";
							}
						print WR "==========  ======  ======  ======  ==========    $tmp_equal_sign    ======\n";
						
						
						my $k=$i+3;
						while($arr_crt[$k]!~/^---/)
							{
								my @arr_line=split('\t+',$arr_crt[$k]);
								
								my $start=$arr_line[0]; 										my $gapped_start=&fill_string_with_gaps($start,10,"LEFT");
								my $repeat_len=length($arr_line[1]);							my $gapped_repeat_len=&fill_string_with_gaps($repeat_len,6,"LEFT");
								my $percent_id="100.0";											my $gapped_percent_id=&fill_string_with_gaps($percent_id,6,"LEFT");
								my $spacer_len=length($arr_line[2]);							my $gapped_spacer_len=&fill_string_with_gaps($spacer_len,6,"LEFT");
								my $left_flank="NNNNNNNNNN";									my $gapped_left_flank=&fill_string_with_gaps($left_flank,10,"LEFT");
								my $repeat=&change_bases_to_dots($arr_line[1],$model_repeat);	my $gapped_repeat=&fill_string_with_gaps($repeat,length($model_repeat),"RIGHT");
								my $spacer="";
								if($arr_crt[$k+1]!~/^---/)
									{
										$spacer=$arr_line[2];
									}
								else{
										$spacer="NNNNNNNNNN";
										$spacer_len=0;
																								$gapped_spacer_len=&fill_string_with_gaps($spacer_len,6,"LEFT");
									}										
																								my $gapped_spacer=&fill_string_with_gaps($spacer,10,"RIGHT");			
								
								
								print WR "$gapped_start  $gapped_repeat_len  $gapped_percent_id  $gapped_spacer_len  $gapped_left_flank    $gapped_repeat    $gapped_spacer\n";
								
								$k++;
							}
						
						#---- now get the last line ----
						my $no_of_repeats;
						my $avg_repeat_len;
						my $avg_p_id=" ";
						my $avg_spacer_len;
						my $avg_left_flank=" ";
						if($arr_crt[$k+1]=~/Repeats:/)
							{
								my @arr_tmp_line=split('\t+',$arr_crt[$k+1]);
								if($arr_tmp_line[0]=~/Repeats: (\d+)/){$no_of_repeats=$1;}else{$no_of_repeats=0;}
								if($arr_tmp_line[1]=~/Length: (\d+)/){$avg_repeat_len=$1;}else{$avg_repeat_len=0;}
								if($arr_tmp_line[2]=~/Length: (\d+)/){$avg_spacer_len=$1;}else{$avg_spacer_len=0;}
							}
						
						my $gapped_no_of_repeats=&fill_string_with_gaps($no_of_repeats,10,"LEFT");	
						my $gapped_avg_repeat_len=&fill_string_with_gaps($avg_repeat_len,6,"LEFT");
						my $gapped_avg_p_id=&fill_string_with_gaps($avg_p_id,6,"LEFT");
						my $gapped_avg_spacer_len=&fill_string_with_gaps($avg_spacer_len,6,"LEFT");
						my $gapped_avg_left_flank=&fill_string_with_gaps($avg_left_flank,10,"LEFT");
						
						print WR "==========  ======  ======  ======  ==========    $tmp_equal_sign\n";						
						print WR "$gapped_no_of_repeats  $gapped_avg_repeat_len  $gapped_avg_p_id  $gapped_avg_spacer_len  $gapped_avg_left_flank    $model_repeat\n";
						print WR "\n";
						print WR "SUMMARY BY SIMILARITY\n";
						close(WR);
						#--- end of current array ----
						$i=$k;
						#exit;
					}
			}
		
		
		
		
		
		
		return 1;
	}




sub fill_string_with_gaps()
	{
		my($object,$length,$side)=@_;
		
		if(not defined $object or $object eq "NA"){$object="";}
		
		my $gap_filled_string=$object;
		
		for(my $i=0;$i<=$length;$i++)
			{
				if(length($gap_filled_string)==$length){last;}
				
				if($side eq "RIGHT")
					{
						$gap_filled_string=$gap_filled_string." ";
					}
				else{
						$gap_filled_string=" ".$gap_filled_string;
					}		
			}
		
		return($gap_filled_string);
	}

sub get_array_degeneracy_score_pcr()
	{
		my($current_repeat,$ref_arr_repeats)=@_;
		
		my $arr_degen_score=0;
		
		foreach my $repeat(@{$ref_arr_repeats})
			{								
				my $sim_score=&get_similarity_score($repeat,$current_repeat);
				
				my $neg_score=length($current_repeat)-$sim_score;
				
				$arr_degen_score=$arr_degen_score-$neg_score;
			}
			
		return($arr_degen_score);	
	}

sub get_similarity_score()
	{
		my($current_repeat,$model_repeat)=@_;
		
		
		
		my $similarity_score=0;
		
		if($current_repeat eq "" or $model_repeat eq ""){return $similarity_score;}
		
		
		my @arr_1=split('',$current_repeat);
		my @arr_2=split('',$model_repeat);
		
		
		
		for(my $i=0;$i<=$#arr_2;$i++)   
			{
				if(not defined $arr_1[$i] or not defined $arr_2[$i]){next;}
				if($arr_1[$i] eq $arr_2[$i])
					{
						$similarity_score=$similarity_score+1;	
					}	
			}
		
		return $similarity_score;
	}


sub change_bases_to_dots()
	{
		my($bf_string,$mr_string)=@_;
		
		if(not $bf_string){return("");}
		#print "\$bf_string,\$mr_string $bf_string,$mr_string\n";
		my $return_string;
		
		my @arr_1=split('',$bf_string);
		my @arr_2=split('',$mr_string);
		
		for(my $i=0;$i<=$#arr_2;$i++)
			{
				if(not defined $arr_1[$i])
					{
						$return_string=$return_string."-";
					}
				elsif($arr_1[$i] eq $arr_2[$i])
					{
						if($arr_2[$i] eq "-")
							{
								$return_string=$return_string."-";
							}
						else{
								$return_string=$return_string.".";
							}
						
					}
				else{
						$return_string=$return_string.$arr_1[$i];
					}	
			}
		
		return($return_string);
	}


	


#-- end ---------------------------------
return 1;


