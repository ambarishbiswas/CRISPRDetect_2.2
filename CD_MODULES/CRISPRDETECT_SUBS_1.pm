#--- begin ----------------------------


################################################ sub routines ###############################################		




sub get_crisprdetect_hotspots()
	{
		my($minimum_repeat_length,$minimum_no_of_repeats,$minimum_spacer_gap,$maximum_spacer_gap,$seq_file,$combined_hotspots_file,$hash_id_lookup_table)=@_;
		
		

		
		
		
		my $id_line=`grep '>' $tmp_dir\/$seq_file >&1`;
		chomp $id_line;$id_line=~s/\r//;$id_line=~s/^>//;
		my $accession=$id_line;
		my $defination=" ";
		if($hash_id_lookup_table->{$accession})
			{
				$defination=$hash_id_lookup_table->{$accession};		
				$defination=~s/^$accession\-?//;
			}	

		
		#print "$id_line\n";
		
		
		
		my $seq=`tail -1 $tmp_dir\/$seq_file >&1`;
		
		#---- get the CRISPR hotspots
		my %hash_of_crispr_hotspots;
		
		#my $method="SLOW";
		my $method="FAST";
		
		if($method eq "FAST")
			{
				while($seq=~/(\w{$minimum_repeat_length})/g)
					{
						my $match_position=$-[0];
						my $seq1=substr($seq,$-[0],1500);
						if($seq1=~/^($1)\w{21,130}?$1/)
						{
							#print "$1\t$match_position\n";
							$hash_of_crispr_hotspots{$match_position}=$1;
						}
					}	
			}
		else{		
				#print "$seq\n";
				my $z=0;
				while($seq=~/\w{$z}(\w{$minimum_repeat_length})/)				
					{
							my $match_position=$z+$-[0];							
							#print "$match_position\t$1\n";
							my $seq1=substr($seq,$-[0],1500);
							if($seq1=~/^\w{$z}($1)\w{21,190}?$1/)							
								{
									#print "$1\t$match_position\n";
									$hash_of_crispr_hotspots{$match_position}=$1;
								}	
						$z=$z+int($minimum_repeat_length/4);
						if($z>=length($seq)){last;}		
					}
			}	
			
		my $total_hotspots=keys %hash_of_crispr_hotspots;
		

		
		
			
			my $seq_len=length($seq);
			my $my_output_file=$accession.&get_unique_id()."_CD_output.txt";

			open(HSP,">$tmp_dir\/$my_output_file") or print "$!";	close(WR);		#system("chmod 777 $tmp_dir\/$my_output_file");
			
			open(HSP,">$tmp_dir\/$my_output_file");
			flock(HSP,2);
			print HSP "ORGANISM:  $accession| \n";
			print HSP "Bases: $seq_len\n";
			#my $str=substr($seq,10,100);
			#print "$str\n";
			my $arr_count=1;
			#my $minimum_repeat_length=11;
			my $max_repeat_length=67;
			
			#for(my $i=1;$i<=length($seq);$i++)
			my $position_checked_already=0;
			foreach my $i(sort{$a<=>$b} keys %hash_of_crispr_hotspots)
			{
				#print "OK1: $i\n";
				
				#if($i<1418032){next;} ##3967392
				#if($i>1419900){next;}
				if($i<$position_checked_already){next;}
				
								
				
				my $seq1=substr($seq,$i,1500);
				#print "OK2: $seq1\n";
				
				
				if($seq1=~/^(\w{$minimum_repeat_length})/)
					{
						#my $seq1=$_;
						my $matching_word=$1;
						#
						#print "$matching_word\t$i $-[0]-$+[0]\n";
					
					my @arr_words;
					my $allowed_one_mismatch=0;
					
					if($allowed_one_mismatch==1)
						{
						
						my @arr_tmp1=split('',$matching_word);
						for(my $j=1;$j<$#arr_tmp1;$j++)
							{
								my $reg_ex="";
								
								for(my $k1=0;$k1<=$#arr_tmp1;$k1++)
									{
										if($j==$k1)
											{
												$reg_ex=$reg_ex.".?";
											}
										else{
												$reg_ex=$reg_ex.$arr_tmp1[$k1];
											}		
									}
								push(@arr_words,$reg_ex);	
							}
						}	
					else{
							push(@arr_words,$matching_word);
						}	
					
					
					foreach my $target_word(@arr_words)
							{
								#print "Checking with $target_word..\n";
								my $match_position;
								my $first_repeat_start;
								my $match_found=0;
								
								
								if($minimum_no_of_repeats<=2)
									{
										if($seq1=~/($matching_word)(\w{$minimum_spacer_gap,$maximum_spacer_gap}?)($target_word)(\w{0,$maximum_spacer_gap}?)/g)
											{
												$match_position=$-[0]; #print "$match_position\n";
												$first_repeat_start=$i+$-[0];
												$match_found=1;
												
												$seq1=substr($seq1,$match_position);
											}
									}
								else{
										if($seq1=~/($matching_word)(\w{$minimum_spacer_gap,$maximum_spacer_gap}?)($target_word)(\w{$minimum_spacer_gap,$maximum_spacer_gap}?)($target_word)/g)
											{
												$match_position=$-[0]; #print "$match_position\n";
												$first_repeat_start=$i+$-[0];
												$match_found=1;
												
												$seq1=substr($seq1,$match_position);
											}	
									}
									
								#if($seq1=~/($matching_word)(\w{25,100})($target_word)(\w{0,100})/g)
								if($match_found==1)
									{
										#print "Match found\n\n";
										#$match_position=$-[0]; #print "$match_position\n";
										#$first_repeat_start=$i+$-[0];

										my %hash_of_potential_crispr;
										my($array_start,$array_stop,$array_length)=&find_initial_crispr($accession,$first_repeat_start,$maximum_spacer_gap,$matching_word,$seq1,\%hash_of_potential_crispr);
										
										#print "$array_start,$array_stop,$array_length\n";
										if($array_stop==0 and $array_length>0)
											{
												$i=$i+$array_length;										
												if($i>$position_checked_already){$position_checked_already=$i;}
												last;
											}
										elsif($array_length>0)
											{
												#451273-452587
												#my $t_seq=substr($seq,451273,100);print "\n$t_seq\n\n";
												
												#--- convert to 1 based coordinate
												$array_start=$array_start+1;
												$array_stop=$array_stop+1;
												
												print HSP "\n\nCRISPR $arr_count   Range: $array_start - $array_stop\n";
												print HSP "POSITION	REPEAT				SPACER\n";
												print HSP "--------	-------------------------------	-----------------------------------------------------------\n";
												my $l1;
												my $l2;
												my $no_of_repeats=0;
												foreach my $row_count(sort{$a<=>$b} keys %hash_of_potential_crispr)
													{
														my($r_start,$r_seq,$s_seq)=split('\t',$hash_of_potential_crispr{$row_count});
														
														$r_start=$r_start+1;  #-- convert to 1 base coordinate system
														
														if(not defined $s_seq){$s_seq="-";}
														else{
																chomp $s_seq;
																
															}
														$s_seq=~s/-//g;
														
														if($s_seq eq ""){$s_seq="-";}
														
														
														if($s_seq!~/\|/)
															{
																$l1=length($r_seq);
																$l2=length($s_seq);
																print HSP "$r_start\t$r_seq\t$s_seq\t[ $l1, $l2 ]\n";	
																#print "$r_start\t$r_seq\t$s_seq\t[ $l1, $l2 ]\n";	
																
															}
														else{
																print HSP "$r_start\t$r_seq\t\n";
																#print "$r_start\t$r_seq\t\n";
																
															}
														$no_of_repeats++;			
													}
											
												
												print HSP "--------	-------------------------------	-----------------------------------------------------------\n";
												print HSP "Repeats: $no_of_repeats	Average Length: $l1		Average Length: $l2\n";
												#print "Repeats: $no_of_repeats	Average Length: $l1		Average Length: $l2\n";
												
												#---- check this later ---										
												$i=$i+$array_length;
												
												if($i>$position_checked_already){$position_checked_already=$i;}
												#print "\n";		
												$arr_count++;
												last;
												#-------------------------									
												
											}

									}
							}
						
					
					
					}
					
				$i++;			
			}
			
			close(HSP);	
			
			#print "Total CRISPR hotspots identified: $arr_count\n";	
			
			#exit;
			
			&simplify_crt_formatted_arrays_and_append($my_output_file,$combined_hotspots_file);
			
			
			unlink("$tmp_dir\/$my_output_file");
			
			return ($total_hotspots,$arr_count);
		}
	


sub find_initial_crispr()
	{
		my ($accession,$start_position,$max_spacer_length,$matching_word,$seq1,$hash_of_potential_crispr)=@_;
		
		#if($matching_word!~/GGGCGGGTGC/){return (0,0,0);}
		my $skip=1;
		
		if($skip==0)
			{
				#print "\nWord: $matching_word at position:$start_position\n";
				#--- find the longest word in the given string ----
				
				my $biggest_word=$matching_word;
				my $found=0;
				for(my $i=38;$i>10;$i--)
					{
						my $j=0;
						while($seq1=~/^\w{$j}(\w{$i})/)
							{
								my $match_position=$start_position+$-[0];
								#print "$match_position\n";
								my $seq2=substr($seq1,$-[0],1500);
								if($seq2=~/($1)\w{21,110}$1/)
									{
										print "$1\t$match_position\n";
										$biggest_word=$1;
										#$hash_of_crispr _hotspots{$match_position}=$1;
										$found=1;
										last;
									}
								$j++;
								if($j>=length($seq1)){last;}		
							}
						#print "\n";	
						
						if($found==1){last;}
					}
					
					
				my @arr_t11=split($biggest_word,$seq1);
				my $left_flank=shift(@arr_t11);
				$start_position=$start_position+length($left_flank);
				foreach my $spacer(@arr_t11)
					{
						#print "$biggest_word-$spacer\n";
					}
				$matching_word=$biggest_word;	
				#print "\n";
			}
		#exit;
		#
		
		my $min_word_len=int(length($matching_word)/2)-1;
		my @arr_spacers=split($matching_word,$seq1);
		
		shift(@arr_spacers); #--- first one will be always blank as the word starts from the first position
		
		#print "@arr_spacers\n";exit;
		
		#----- check if no_of_clusters =1; if so return -------
		#----- don't include the last spacer either, as they are mostly non-identical or long
		my @tmp_arr_spacers;#=@arr_spacers;
		my %tmp_hash_of_spacers;
		for(my $i=0;$i<=$#arr_spacers;$i++)
			{
				if(length($arr_spacers[$i])<=125)
					{
						$tmp_arr_spacers[$i]=$arr_spacers[$i]; #---- remove the last spacer
						
						if(defined $tmp_hash_of_spacers{$arr_spacers[$i]})
							{
								$tmp_hash_of_spacers{$arr_spacers[$i]}=$tmp_hash_of_spacers{$arr_spacers[$i]}+1;
							}
						else{
								$tmp_hash_of_spacers{$arr_spacers[$i]}=1;
							}	
					}
				else{
						last;
					}		
			}
				
		if(($#tmp_arr_spacers+1) >= 2)
			{
				
				#print "Processing...\n";
				#foreach (@tmp_arr_spacers)
				#	{
				#		print "$_\n";
				#	}
				
				#--------- first check identical spacers --------------------------------------------------------------	
				my $max_identical_spacers=0;	
				foreach my $spacer_t(sort{$tmp_hash_of_spacers{$b}<=>$tmp_hash_of_spacers{$a}} keys %tmp_hash_of_spacers)
					{
						$max_identical_spacers=$tmp_hash_of_spacers{$spacer_t};
						#print "$spacer_t\t$max_identical_spacers\n";
						last;
					}
				if($max_identical_spacers >=(($#tmp_arr_spacers+1)/2+1)) #-- more than 50% spacers are identical
					{
						my $len_t=length($matching_word)+length($tmp_arr_spacers[0])+length($matching_word);
						return($start_position,0,$len_t);
					}
				
					
				#---------- then check no. of clusters -----------------------------------------------------------------
					
				my $total_number_of_clusters=int(&get_spacers_identity($start_position,$accession,\@tmp_arr_spacers));
				#print "\$max_identical_spacers=$max_identical_spacers\t\$total_number_of_clusters=$total_number_of_clusters\n";
				if($total_number_of_clusters<=1)
					{
						my $len_t=length($matching_word)+length($tmp_arr_spacers[0])+length($matching_word);
						return($start_position,0,$len_t);
					}
			}
		
		
		
		
		
		
		#---now check for tandem repeats by joining the word and the first spacer, and further split to check if the region contains long tandem repeats -----
		
		my $skip_2=1;
		
		if($skip_2==0)
		{
		if(defined $arr_spacers[0])
			{
				
				my $checking_word=$matching_word.$arr_spacers[0];
				
				my @arr_spacers2=split($checking_word,$seq1);
				shift(@arr_spacers2);
				my $valid_spacers=1;
				foreach my $spacer2(@arr_spacers2)
					{
						if(length($spacer2)<100){$valid_spacers++}else{last;}
						#print "$checking_word-$spacer2\n";
						
					}
				if($valid_spacers>1)
					{
						my $array_stop=0;
						my $array_length=length($checking_word)*$valid_spacers;
						
						return($start_position,$array_stop,$array_length);
					}
			}
		}
		
		#-----------------------------------------------------------------------------------------------------------------------------------------------------	
			
	
		my $time_1 = "CHP_".$accession.$start_position.&get_unique_id().$matching_word;	   

		my $tmp_file=		$time_1."_tmp_spacers.txt";
		my $tmp_output=		$time_1."_tmp_spacers.aln";
		my $tmp_output_dnd=	$time_1."_tmp_spacers.dnd";
		
		open(FA,">$tmp_dir\/$tmp_file") or print "$!";
		flock(FA,2);
		my $spacer_index=0;
		my $longest_spacer=0;
		my $stop_flag=0;
		foreach my $spacer(@arr_spacers)
			{
				if($spacer_index<$#arr_spacers and length($spacer)>$longest_spacer){$longest_spacer=length($spacer);}
								
				if($spacer_index==$#arr_spacers){$spacer=substr($spacer,0,$longest_spacer);$stop_flag=1;}
				
				if(length($spacer)>$max_spacer_length){$spacer=substr($spacer,0,$max_spacer_length);$stop_flag=1;}
				
				print FA ">S_$spacer_index\n$matching_word$spacer\n";
				
				if($stop_flag==1){last;}
				$spacer_index++;
			}
		close(FA);
		
		#---- run clustalW to get the repeat and spacers ------------------------------------------------------------------------------------------------
		
		my @arr_output;
		my %hash_of_lines;
		
		#print "clustalw -INFILE=$tmp_dir\/$tmp_file -OUTFILE=$tmp_dir\/$tmp_output -QUIET -ALIGN -GAPOPEN=0.01 -GAPEXT=0.01 -OUTORDER=INPUT\n";
		system("clustalw -INFILE=$tmp_dir\/$tmp_file -OUTFILE=$tmp_dir\/$tmp_output -QUIET -ALIGN -NUMITER=50 -NOWEIGHTS -ENDGAPS -OUTORDER=INPUT >/dev/null 2>&1");
		
		if(-e "$tmp_dir\/$tmp_output")
			{
				open(RD,"$tmp_dir\/$tmp_output");
				@arr_output=<RD>;
				close(RD);
				
				unlink("$tmp_dir\/$tmp_file");
				unlink("$tmp_dir\/$tmp_output");
				unlink("$tmp_dir\/$tmp_output_dnd");
			}
		#--------------------------------------------------
		
		
		for(my $i=0;$i<=$#arr_output;$i++)
			{
				my $line=$arr_output[$i];
				chomp $line; $line=~s/\r$//;
				#print "$line\n";
				if($line=~/^S_/)
					{	
											
						
						$line=~s/^S_//; $line=~s/\s+/\t/g;
						my @arr_t1=split('\t',$line);
						
						if(defined $hash_of_lines{$arr_t1[0]})
							{
								$hash_of_lines{$arr_t1[0]}=$hash_of_lines{$arr_t1[0]}.$arr_t1[1];
							}
						else{
								
								$hash_of_lines{$arr_t1[0]}=$arr_t1[1];
							}	
					}
				elsif($line=~/^\s{16}/)
					{
						$line=~s/^\s{16}//;
						
						my $hash_index=$#arr_spacers+1;
						if($arr_output[$i-1]=~/^S_(\d+)/){$hash_index=$1+1;}
						
						if(defined $hash_of_lines{$hash_index})
							{
								$hash_of_lines{$hash_index}=$hash_of_lines{$hash_index}.$line;
							}
						else{								
								if($line ne "")
									{
										$hash_of_lines{$hash_index}=$line;
									}	
							}
					}	
					
			}
			
		#my $op=system("cat $tmp_dir\/$tmp_output >&1");
		
		#print "$op\n";
		my $exit_code=0;
		my $starting_gap=0;
		my $model_repeat="";
		foreach my $seq_index(sort{$a<=>$b} keys %hash_of_lines)
			{
				#print "$seq_index\t$hash_of_lines{$seq_index}\n";
				
				if($hash_of_lines{$seq_index}=~/\*/)
					{
						if($hash_of_lines{$seq_index}=~/^(\s{0,$max_spacer_length})(\*{$min_word_len,$max_spacer_length})\s{0,}/)
							{
								if(defined $1){$starting_gap=length($1);}
								
								#print "R $#arr_spacers:\t$2\n\n";
								$model_repeat=$2; $model_repeat=~s/\s+$//;								
							}
						
						#$exit_code=1;
					}
			}
		if(not defined $model_repeat or length($model_repeat)<length($matching_word)){return(0,0,0);}
		#----- now add the bases falling in the initial gap
		my $aditional_bases_at_the_beginning="";
		if($starting_gap>0)
			{
				foreach my $seq_index(sort{$a<=>$b} keys %hash_of_lines)
					{
						#print "$seq_index\t$hash_of_lines{$seq_index}\n";
						my $line=$hash_of_lines{$seq_index};
						
						my $leading_bases=substr($line,0,$starting_gap);
						if($seq_index==0)
							{
								$aditional_bases_at_the_beginning=$leading_bases; $aditional_bases_at_the_beginning=~s/-//g;
							}
						else{
								#if(not defined $hash_of_lines{$seq_index-1}){next;}
								$hash_of_lines{$seq_index-1}=$hash_of_lines{$seq_index-1}.$leading_bases;
							}	
						
						my $current_line=substr($line,$starting_gap);
						$hash_of_lines{$seq_index}=$current_line;
					}
					
							
			}
			
		$start_position=$start_position+length($aditional_bases_at_the_beginning);
		
		#---- now split the repeat and spacer for each record
		my $repeat_start=0;
		my $array_stop=0;
		my $array_length=0;
		my $total_spacer_length=0;
		my $number_of_repeats=0;
		foreach my $seq_index(sort{$a<=>$b} keys %hash_of_lines)
			{
				if(not defined $hash_of_lines{$seq_index} or $hash_of_lines{$seq_index}=~/\*/){next;}
				my $line=$hash_of_lines{$seq_index};
				
				my $repeat=substr($line,0,length($model_repeat)); my $repeat_1=$repeat; $repeat_1=~s/-//g;
				my $spacer=substr($line,length($model_repeat));$spacer=~s/-//g;
				
				$total_spacer_length=$total_spacer_length+length($spacer);
				
				if($seq_index==0){$repeat_start=$start_position;}
				else{
						#if(not defined $hash_of_lines{$seq_index-1}){next;}
						
						my $previous_line=$hash_of_lines{$seq_index-1};$previous_line=~s/-//g;
						
						$repeat_start=$repeat_start+length($previous_line);
					}
				
				#print "$repeat_start\t$repeat\t$spacer\n";	
				if($hash_of_lines{$seq_index+1}=~/\*/)
					{
						$spacer="|";
						#print "$repeat_start\t$repeat\t|\n";
						$array_stop=$repeat_start+length($repeat);
						$array_length=$array_length+length($repeat_1);
					}
				else{
						$line=~s/-//g;
						$array_length=$array_length+length($line);
					}	
				
				
				$hash_of_potential_crispr->{$seq_index}="$repeat_start\t$repeat\t$spacer";		
				
				$number_of_repeats++;		
			}
		
		my $avg_spacer_len=int($total_spacer_length/($number_of_repeats-1));
		#$array_length=$array_stop-$start_position+1;
		
		unlink<$tmp_dir\/$time_1\*>;
		
		
		if($avg_spacer_len<5){return($start_position,0,0);}
		#if(length($model_repeat)<17 and $avg_spacer_len>70){return($start_position,0,0);} #--- block it later. bad way to suppress possible CRISPRS
		
		return($start_position,$array_stop,$array_length);
			
	}



sub load_crispr_hotspots()
	{		
		my($input_array_file,$hash_of_all_crispr_hotspots)=@_;
		
		

		open(RD,$input_array_file) or print "$!";
		my @arr_compiled_pilercr_output=<RD>;
		close(RD);

		my @array_cpr_lines=@arr_compiled_pilercr_output;


		#my $cnt_2=1;
		my $cnt_2_2=1;
		my $cnt_2_3=1;



		my $array_index=0;	
		my $total_records=0;



		#while(my $line=<RD>)
		for(my $i=0;$i<=$#arr_compiled_pilercr_output;$i++)
			{		
				my $line=$arr_compiled_pilercr_output[$i];
				
				if($line=~/^>/ and $arr_compiled_pilercr_output[$i+2]=~/Pos/)
					{
						#print "$line\n<br>";
						my($accession,$repeat_sequence2,$no_of_spacers2);
						my($pcr_crispr_start,$pcr_crispr_stop,$repeat_start,$repeat_stop);
						
				
						my @tmp_arr_id=split('\|',$line);
						my @current_array;	
						
						my $array_line_start=$i;
								
						$accession=$tmp_arr_id[0]; chomp $accession;$accession=~ s/^>//; 
						
									#-- this is for cases where there are two accessions, one old and one new: 
													#DEFINITION  Candidatus Cloacamonas acidaminovorans str. Evry provisional genome
													#			sequence from WWE1 candidate division.
													#ACCESSION   NC_020449 NS_000195
						if($accession=~/^(\S+)/)  
							{
								$accession=$1;
							}
						$accession=~ s/\.\d+$//; $accession=~s/\r//;
						

						
						#################################################################
						
						
						$array_index++;

						
						if($arr_compiled_pilercr_output[$i+2]!~/\s+Pos\s+Repeat/){next;}
							
						
						
						my $coord_line=$arr_compiled_pilercr_output[$i+4]; 	#   1372383      26    96.2      43  ACTGGAGAGT    ......C...................    CAATTTAGAATGGCTACAAGCCGATGGTAATCAGCTAAGTCGG
						
						$coord_line=~s/^\s+//;
						
						my @arr_tmp=split('\s+',$coord_line);
						
						$pcr_crispr_start=$arr_tmp[0];
						$repeat_start=$arr_tmp[0];
						
					
						#------- get the last repeat -------------------------------------------------------------------------------------------
						
						my $j=$i+4;
						while($arr_compiled_pilercr_output[$j]!~/=======/)
							{
								if($arr_compiled_pilercr_output[$j+1]=~/^====/)#  10095660      36   100.0          TGGCAAACCA    ....................................    GTACTT
									{
										my $last_line=$arr_compiled_pilercr_output[$j];
										$last_line=~s/^\s+//;$last_line=~s/\s+/\t/g;
										if($last_line=~/^(\d+)\t(\d+)\t/)
											{
												$pcr_crispr_stop=$1+$2-1;
											}
									}	
			
								$j++;	
							}
				

				
				#-------- now get the model repeat seq --------------------------
				my $pcr_crispr_range=$pcr_crispr_start."-".$pcr_crispr_stop;
				my $range="$pcr_crispr_start-$pcr_crispr_stop";
						
				my $array_line_stop=$j+1;		
							
					
				

						
				#----- now store the whole array in @original_array		
				
				
				my @original_pcr_array;
				
				for(my $k=$array_line_start;$k<=$array_line_stop;$k++)
					{
						#print "$arr_compiled_pilercr_output[$k]";
						chomp $array_cpr_lines[$k]; $array_cpr_lines[$k]=~ s/\r+//g;
						push(@original_pcr_array,$array_cpr_lines[$k]);
					}	
				
				my $orig_array=join("\n",@original_pcr_array);
					
				$hash_of_all_crispr_hotspots->{$accession}->{$range}=$orig_array;
				
			
				#print "$i - $array_line_stop: \$hash_of_all_crispr_hotspots{$accession}{$range}\n$hash_of_all_crispr_hotspots{$accession}{$range}\n\n";
				$i=$array_line_stop;
				#next;
				
				
				
			}
		 else{next;}
			
		 #	if($single_array_test>20){last;}	
		}		

		return($total_records);
	}		
	



sub extract_information_and_sequence_from_gbk_file()
	{
		my($filename)=shift(@_);
		
		
		my $accession=`grep -w "ACCESSION  " $filename >&1`;  #---- do not change grep condition blindly, the No. of gapa are important
		$accession=~s/ACCESSION  //;chomp $accession; $accession=~s/\r//g;$accession=~s/^\s+//;$accession=~s/\s+$//;
		
		if($accession=~/\s+/)
			{
				#---- there are two or more accessions like "NC_014820 NS_000189"
				my @arr_t1=split(' ',$accession);
				$accession=$arr_t1[0];
			}		
				
				
		my $defination=`grep -w "DEFINITION " $filename >&1`;
		$defination=~s/DEFINITION //; chomp $defination;$defination=~s/\r//g;$defination=~s/^\s+//;
				
				#print "\$accession=$accession\n\$defination=$defination\n";exit;
				
		my $species="Bacteria";		
		my $species_line=`cat $filename | grep -A 2 -w 'ORGANISM' | grep -E 'Bacteria|Archaea|Viruses' >&1`;
		if(defined $species_line and $species_line ne "")
			{
				chomp $species_line;$species_line=~s/\r//g;$species_line=~s/^\s+//;
				my @arr_species_line=split(';',$species_line);
						
						#print "$arr_species_line[0]\n\n";
						
				$species=$arr_species_line[0];		
				#$hash_of_accession_and_species{$accession}=$arr_species_line[0];
			}
		#----------- call seqret to extract the fasta sequence-----------------------------------------
				
		#if(not -e "$tmp_dir\/$accession\.fna" or -s "$tmp_dir\/$accession\.fna" < 10)
		#	{												
						my $tmp_file=$accession.&get_unique_id()."_tmp_seq.txt";
						
						#print "";
						
						if($filename=~/\:/) #----- this causes problem with seqret; create a .gbk file and replace the $filename to fool seqret
							{
								my $tmp_f1="$tmp_dir\/".$accession.".gbk";
								system("cp $filename $tmp_f1");
								$filename=$tmp_f1;
							}
						
						my $ret_text=`seqret -sequence $filename -outseq $tmp_dir\/$tmp_file >/dev/null 2>&1`;
						#my $ret_text=`perl gbk2fasta.pl $filename $tmp_dir\/$tmp_file >/dev/null 2>&1`;
						
						my $file_size=-s "$tmp_dir\/$tmp_file";
						if(not defined $file_size or $file_size<=0)
							{	
								my $ret_text2=`seqret -sequence $filename -outseq $tmp_dir\/$tmp_file >/dev/null 2>&1`;
									
								my $file_size2=-s "$tmp_dir\/$tmp_file";							
								if(not defined $file_size2 or $file_size2<=0)
									{
										print "Error: File size=$file_size2 The GBK file $filename does not contain any genomic sequence. Please download a full .gbk file\n"; return('NA',$defination,$species);
									}	
							}
						#---- now open the $tmp_file and convert it to supported .fna style
						my $seq="";
						open(RD,"$tmp_dir\/$tmp_file") or print "$!";
						flock(RD,2);
						while(my $line=<RD>)
							{
								if($line!~/>/)
									{
										chomp $line;
										$line=~s/\r//g;
										$line=uc($line);
										$seq=$seq.$line;
									}
							}
						close(RD);
						
						#----- now write the sequence to $accession.fna file
						open(WR,">$tmp_dir\/$accession\.fna") or print "$!";
						flock(WR,2);
						print WR ">$accession\n$seq\n";
						close(WR);
						
						if(-e "$tmp_dir\/$accession\.fna"){system("chmod 777 $tmp_dir\/$accession\.fna");}
						unlink("$tmp_dir\/$tmp_file");	
			#}
			
		return($accession,$defination,$species);	
			
				
	}



sub split_gbff_file_to_individual_gbk_files()
	{
		my($input_gbff_file,$arr_gbk_files)=@_;
		
		#print "\n$input_gbff_file\n";
		my @arr_line=`less $input_gbff_file >&1`;
		
		my $last_acc="";
		for(my $i=0;$i<=$#arr_line;$i++)
			{
				#print "$arr_line[$i]";exit;
				if($arr_line[$i]=~/^LOCUS    /)
					{
						my @arr_t1=split('\s+',$arr_line[$i]);
						$last_acc=$arr_t1[1];
						
						
						my $j=$i;
						my $gbk_file="$tmp_dir\/".$last_acc.".gbk";
						open(WR,">$gbk_file");
						
						while($arr_line[$j]!~/^\/\//)
							{
								print WR "$arr_line[$j]";
								
								
								$j++;
								if($j>=$#arr_line){last;}
							}
							
						print WR "$arr_line[$j]";	
						close(WR);
						push(@{$arr_gbk_files},$gbk_file);
						#print "\t\$gbk_file= $gbk_file\n";
						$i=	$j;
						
					}
			}
		
		return 1;
		
	}
sub process_gbk_file()
	{
		my($input_gbk_file,$output_gene_position_file)=@_;
		
		#my @arr_tmp_1=split("/",$input_gbk_file);		
		#my $accession=$arr_tmp_1[$#arr_tmp_1]; $accession=~s/\.\S+$//;
		
		#print "Accession : $accession [$input_gbk_file] \$input_gbk_file=$input_gbk_file \n";

		#---- first load Cas gene synonyms ------------------------
		my @arr_cas_synonym_lines=`less Ref_lib_files/cas_gene_synonyms.tab | grep -v '^#' >&1`;
		my %hash_of_synonyms_and_proposed_cas_name;

		foreach my $line(@arr_cas_synonym_lines)
			{
				chomp $line;$line=~s/\r//g;
				if($line!~/\S+/){next;}
				
				my @arr_t1=split('\t',$line);
				my $proposed_cas_name=$arr_t1[0];
				my $synonyms=$arr_t1[1];
				my @arr_synonyms=split(',',$synonyms);
				$hash_of_synonyms_and_proposed_cas_name{$proposed_cas_name}=$proposed_cas_name;
				
				foreach my $s_name(@arr_synonyms)
					{
						$hash_of_synonyms_and_proposed_cas_name{$s_name}=$proposed_cas_name;
					}
			}

		
		#-------------------------------------------------------------------------------
		open(RD,"$input_gbk_file") or print "$!";		
		my @arr_gbk_file=<RD>;
		close(RD);
		
		open(TAB,">>$tmp_dir\/$output_gene_position_file") or print "";
		flock(TAB,2);
		
		my @arr_info_line=`grep -A 2 -E '^LOCUS|^VERSION|strain\=|^     CDS|CRISPR|protein_id\=' $input_gbk_file  | grep -v '^-' >&1`;
		
		my $current_accession="";
		my $cas_gene_start;
		my $cas_gene_stop;
		
		my %hash_of_accession_and_gi;
		my %hash_of_gi_and_accession;
		for(my $i=0;$i<=$#arr_info_line;$i++)
			{
				#chomp $line; $line=~s/\r//;
				
				
				#----- get the accession and strain [if present]
				if($arr_info_line[$i]=~/^LOCUS    /)
					{
						my @arr_t1=split('\s+',$arr_info_line[$i]);
						$current_accession=$arr_t1[1];
					}	
						
						
						#---- get the version info -----------------------------------------
						
						
				elsif(defined $arr_info_line[$i] and $arr_info_line[$i]=~/VERSION/)
					{
								my $version_acc="NA";
								my $version_gi="NA";
								my $v_line=$arr_info_line[$i];
								
								my @arr_t2=split('\s+',$v_line);
								$version_acc=$arr_t2[$#arr_t2-1];chomp $version_acc;$version_acc=~s/\r//;
								$version_gi=$arr_t2[$#arr_t2];chomp $version_gi;$version_gi=~s/\r//; $version_gi=~s/GI://;
								
								print TAB "$current_accession\tVERSION\t0\t0\t$version_acc\t$version_gi\n";	
								$hash_of_accession_and_gi{$current_accession}=$version_gi;
								$hash_of_gi_and_accession{$version_gi}=$current_accession;
					}
						
							
						#---- get the strain info ------------------------------------------
						
				elsif(defined $arr_info_line[$i] and $arr_info_line[$i]=~/strain\=/)
					{
								my $strain="NA";
								my @arr_t2=split('strain\=',$arr_info_line[$i]);
								$strain=$arr_t2[$#arr_t2];chomp $strain;$strain=~s/\r//; $strain=~s/\"//g;
								print TAB "$current_accession\tSTRAIN\t0\t0\t$strain\t$strain\n";
					}
				#------ now get CDS and CRISPRs
				elsif($arr_info_line[$i]=~/^     CDS/)
					{
						$arr_info_line[$i]=~s/>//g;$arr_info_line[$i]=~s/<//g;
						
						
						if($arr_info_line[$i]=~/(\d+\.\.\d+)/)
							{
								my($start,$stop)=split('\.\.',$1);
								$cas_gene_start=$start;
								$cas_gene_stop=$stop;
								print TAB "$current_accession\tCDS\t$cas_gene_start\t$cas_gene_stop\tNA\tNA\n";
							}
					}	
				elsif(defined $arr_info_line[$i] and $arr_info_line[$i]=~/CRISPR/ and $arr_info_line[$i]=~/\/product\=/)
					{
								
								my $rec_line=$arr_info_line[$i]; chomp $rec_line; $rec_line=~s/\r//;
								
								#---- check if the rec_line is complete or not [complete lines will end with \"] 
								if($rec_line!~/\"$/)
									{
										my $nxt_line=$arr_info_line[$i+1];
										chomp $nxt_line; $nxt_line=~s/\r//g; $nxt_line=~s/^\s+//;
										if($nxt_line!~/\=/)
											{
												$rec_line=$rec_line." ".$nxt_line;
												#print "Joined: $rec_line\n";
											}
									}
								
								
								#--------- get the Cas protein or nuclease type ---------------------	
								my $cas_gene_family_found=0;
								my $cas_det_line="";	
								my $cas_family="Unclassified_Cas_protein";
								my $cas_protein_acc="NA";
								
								my @arr_tmp1=split('\"',$rec_line);
								
								my $info="Unclassified_Cas_protein";
								
								if(defined $arr_tmp1[1] and $arr_tmp1[1]=~/\S+/)
									{
										$info=$arr_tmp1[1];
									}
								elsif(defined $arr_tmp1[0] and $arr_tmp1[0]=~/\S+/)
									{
										$info=$arr_tmp1[0];
										#print "$rec_line\n";
									}	
								else{
										#print "Error: couldn't locate details from $rec_line\n";
									}
								
								
								
								
								#--- remove some unwanted extra info example # CRISPR-associated protein Cas10/Csm1, subtype III-A/MTUBE; Region: cas_TM1811_Csm1; TIGR02578
								if($info=~/; /)
									{
										#$info=~s/$1//;
										my @arr_tmp1=split(';',$info);
										$info=$arr_tmp1[0];
										#print "Removing: $1\t from: $rec_line\n";
										#print "\t\t$info\n";
									}
								if($info=~/(subtype \w+\-\w+)/)
									{
										$info=~s/$1//;
									}
								
								if($info=~/(Hmari subtype)/)
									{
										$info=~s/$1//;
									}	
								if($info=~/(MYXAN subtype)/)
									{
										$info=~s/$1//;
									}
								if($info=~/(Cyano-type)/)
									{
										$info=~s/$1//;
									}	
								if($info=~/(Anaes-subtype)/)
									{
										$info=~s/$1//;
									}
								if($info=~/(Tneap subtype)/)
									{
										$info=~s/$1//;
									}	
								if($info=~/(subtypeI-B\/HMARI)/)
									{
										$info=~s/$1//;
									}
								if($info=~/(type TIGR02547-like protein)/)
									{
										$info=~s/$1//;
									}					
								#---- remove type and subtype from info
								if($info=~/ (subtype \S+)$/)
									{
										$info=~s/$1//;
										#print "Removing: $1\t from: $rec_line\n";
										#print "\t\t$info\n";
									}
								elsif($info=~/ (type \S+)$/)
									{										
										$info=~s/$1//;
										#print "Removing: $1\t from: $rec_line\n";
										#print "\t\t$info\n";
									}	
								
								if($info=~/,(\S+\s?subtype)$/)
									{
										$info=~s/$1//;
									}
									
								if($info=~/(subtype TIGR01907-like protein)/)
									{
										$info=~s/$1//;
									}	
								
								if($info=~/type/ or $info=~/subtype/)
									{
										#print "\tNeed to process this: $info\n\n";
									}	
									
									
								$info=~s/\(cas\d?\)//gi;		
								$info=~s/\(//g;
								$info=~s/\)//g; # do it twice as some time nested
										
								$info=~s/plasmid//gi;
										
										
								$info=~s/CRISPR-associated//gi;
								$info=~s/CRISPR- associated//gi;
								
								$info=~s/CRISPR system//gi;
								
								$info=~s/family//gi;
								$info=~s/protein//gi;	
								$info=~s/MULTISPECIES://gi;	
								$info=~s/multi-domain//gi;							
								$info=~s/domain//gi;
								
										
								$info=~s/subtype//gi;
								$info=~s/type//gi;
								$info=~s/Cas_//gi;
								$info=~s/Cas\///gi;
									
								$info=~s/^\s+//;
								$info=~s/\s+$//;		
								$info=~s/^,//g;
								$info=~s/,$//g;
								$info=~s/\/\s+/\//g;
								$info=~s/\s+/ /g;
								
								
								#print "$info\t$rec_line\n";	
								
								
								#--- now replace the info with proposed name and add the other info next to the proposed name
								#--- check if multiple proposed name is found suggesting fused Cas genes
								my %hash_of_proposed_names;
								my $proposed_name_found=0;
								foreach my $s_name(sort{length($b)<=>length($a)}keys %hash_of_synonyms_and_proposed_cas_name)
									{
										if($info=~/\b($s_name)\b/i)
											{
												my $m_name=$1;
												my $proposed_cas_name=$hash_of_synonyms_and_proposed_cas_name{$s_name};												
												$hash_of_proposed_names{$proposed_cas_name}=1;
												
												#---- replace the cas name with proposed cas_name
												$info=~s/$m_name/$proposed_cas_name/g;
												
												$proposed_name_found++;
											}
									}
								my $nof_proposed_cas_genes_found=keys %hash_of_proposed_names;
								
								#------------------------------- final cleaning up ---------
								
								$info=~s/-containing//;
								$info=~s/homolog//;
								$info=~s/^\s+//;
								$info=~s/\s+$//;
								
								$info=~s/\s+,\s+/,/;
								$info=~s/\s+\/\s+/\//;
								$info=~s/\s+\//\//;
								$info=~s/\/\s+/\//;
								
								$info=~s/-$//;
								$info=~s/^-//;
								
								$info=~s/\/$//;
								$info=~s/^\///;
								
								if($info=~/(\S+)\/(\S+)\/(\S+)/)
									{
										my $w1=$1;
										my $w2=$2;
										my $w3=$3;
										if($w1 eq $w2 and $w2 eq $w3)
											{
												$info=~s/$w1\/$w2\/$w3/$w1/;
											}	
									}
								if($info=~/(\S+)\/(\S+)/)
									{
										my $w1=$1;
										my $w2=$2;
										if($w1 eq $w2)
											{
												$info=~s/$w1\/$w2/$w1/;
											}	
									}									
									
								$info=~s/,/ /g;	
								$info=~s/-/ /g;	
								$info=~s/\s+/ /g;
								$info=~s/ramp/RAMP/g;
								$info=~s/SAG0894/Cas9/g;
								$info=~s/Cas3\'/Cas3/gi;
								#-----------------------------------------------------------
								if($info eq "Cas" or $info eq "cas")
									{
										$info="";
									}
								if(defined $info and $info=~/\S+/)
									{
										if($proposed_name_found>0)
											{
												$cas_family=$info;
											}
										else{
												$cas_family=$info." Unclassified_Cas_protein";
											}		
									}
								#else{
								#		
								#	}	
								
																
								if($nof_proposed_cas_genes_found>1)
									{
										#print "Multiple proposed_cas_gene_found in $rec_line\n\t$cas_family\n\n";
									}	
								#------	
									
								my $skip_t1=1;
								if($skip_t1==0)
								{	
								my @arr_t2=split('CRISPR-associated ',$rec_line);
								my $cas_info=$arr_t2[1]; $cas_info=~s/\"//g;
								my @arr_t3=split('\s+',$cas_info);	
								if(defined $arr_t3[1] and defined $arr_t3[1]=~/\S+/)
									{
										$arr_t3[1]=~s/,//g;
										$arr_t3[1]=~s/;//g;
										
										$arr_t3[1]=~s/-/_/g;
										
										if($arr_t3[0]=~/protein/i)
											{
												if($arr_t3[1]=~/\S+/)
													{
														$cas_family=$arr_t3[1]."-".$arr_t3[0];;
													}
												else{
														$cas_family="Unclassified_Cas_protein";
													}	
											}
										elsif($arr_t3[1]=~/family/i)
											{
												$cas_family=$arr_t3[0];
											}							
										else{
												$cas_family=$arr_t3[1]."-".$arr_t3[0];
											}	
									}	
								}	
								#----- get the cas_protein_accession
								if(defined $arr_info_line[$i+1] and $arr_info_line[$i+1]=~/protein_id\=/)
									{
										$cas_protein_acc=$arr_info_line[$i+1];chomp $cas_protein_acc; $cas_protein_acc=~s/\r//g; 
										$cas_protein_acc=~s/\/protein_id\=//; $cas_protein_acc=~s/^\s+//; $cas_protein_acc=~s/\"//g; $cas_protein_acc=~s/\.\d+$//;
									}
								elsif(defined $arr_info_line[$i+2] and $arr_info_line[$i+2]=~/protein_id\=/)
									{
										$cas_protein_acc=$arr_info_line[$i+2];chomp $cas_protein_acc; $cas_protein_acc=~s/\r//g; 
										$cas_protein_acc=~s/\/protein_id\=//; $cas_protein_acc=~s/^\s+//; $cas_protein_acc=~s/\"//g; $cas_protein_acc=~s/\.\d+$//;
									}	
																	
								#if($cas_protein_acc eq "NA")
								#	{
								#		print "\n\nError processing: $input_gbk_file\n\t$current_accession\tCRISPR\t$cas_gene_start\t$cas_gene_stop\t$cas_family\t$cas_protein_acc\n";exit;
								#	}	
								print TAB "$current_accession\tCRISPR\t$cas_gene_start\t$cas_gene_stop\t$cas_family\t$cas_protein_acc\n";
						}			
				
				#--------- for older gbk files ------------------
				elsif(defined $arr_info_line[$i] and $arr_info_line[$i]=~/CRISPR/ and $arr_info_line[$i]=~/\/note\=\"CRISPR/)
					{
								
								my $rec_line=$arr_info_line[$i]; chomp $rec_line; $rec_line=~s/\r//;
								if($rec_line=~/CRISPR repeat /){next;}
								#---- check if the rec_line is complete or not [complete lines will end with \"] 
								if($rec_line!~/\"$/)
									{
										my $nxt_line=$arr_info_line[$i+1];
										chomp $nxt_line; $nxt_line=~s/\r//g; $nxt_line=~s/^\s+//;
										if($nxt_line!~/\=/)
											{
												$rec_line=$rec_line." ".$nxt_line;
												#print "Joined: $rec_line\n";
											}
									}
								
								
								#--------- get the Cas protein or nuclease type ---------------------	
								my $cas_gene_family_found=0;
								my $cas_det_line="";	
								my $cas_family="Unclassified_Cas_protein";
								my $cas_protein_acc="NA";
								
								my @arr_tmp1=split('\"',$rec_line);
								
								my $info="Unclassified_Cas_protein";
								
								if(defined $arr_tmp1[1] and $arr_tmp1[1]=~/\S+/)
									{
										$info=$arr_tmp1[1];
									}
								elsif(defined $arr_tmp1[0] and $arr_tmp1[0]=~/\S+/)
									{
										$info=$arr_tmp1[0];
										#print "$rec_line\n";
									}	
								else{
										#print "Error: couldn't locate details from $rec_line\n";
									}
								
								
								
								
								#--- remove some unwanted extra info example # CRISPR-associated protein Cas10/Csm1, subtype III-A/MTUBE; Region: cas_TM1811_Csm1; TIGR02578
								if($info=~/; /)
									{
										#$info=~s/$1//;
										my @arr_tmp1=split(';',$info);
										$info=$arr_tmp1[0];
										#print "Removing: $1\t from: $rec_line\n";
										#print "\t\t$info\n";
									}
								if($info=~/(subtype \w+\-\w+)/)
									{
										$info=~s/$1//;
									}
								
								if($info=~/(Hmari subtype)/)
									{
										$info=~s/$1//;
									}	
								if($info=~/(MYXAN subtype)/)
									{
										$info=~s/$1//;
									}
								if($info=~/(Cyano-type)/)
									{
										$info=~s/$1//;
									}	
								if($info=~/(Anaes-subtype)/)
									{
										$info=~s/$1//;
									}
								if($info=~/(Tneap subtype)/)
									{
										$info=~s/$1//;
									}	
								if($info=~/(subtypeI-B\/HMARI)/)
									{
										$info=~s/$1//;
									}
								if($info=~/(type TIGR02547-like protein)/)
									{
										$info=~s/$1//;
									}					
								#---- remove type and subtype from info
								if($info=~/ (subtype \S+)$/)
									{
										$info=~s/$1//;
										#print "Removing: $1\t from: $rec_line\n";
										#print "\t\t$info\n";
									}
								elsif($info=~/ (type \S+)$/)
									{										
										$info=~s/$1//;
										#print "Removing: $1\t from: $rec_line\n";
										#print "\t\t$info\n";
									}	
								
								if($info=~/,(\S+\s?subtype)$/)
									{
										$info=~s/$1//;
									}
									
								if($info=~/(subtype TIGR01907-like protein)/)
									{
										$info=~s/$1//;
									}	
								
								if($info=~/type/ or $info=~/subtype/)
									{
										#print "\tNeed to process this: $info\n\n";
									}	
									
									
								$info=~s/\(cas\d?\)//gi;		
								$info=~s/\(//g;
								$info=~s/\)//g; # do it twice as some time nested
										
								$info=~s/plasmid//gi;
										
										
								$info=~s/CRISPR-associated//gi;
								$info=~s/CRISPR- associated//gi;
								
								$info=~s/CRISPR system//gi;
								$info=~s/CRISPR\/Cas system-associated//gi;
								
								$info=~s/superfamily//gi;
								$info=~s/family//gi;
								$info=~s/protein//gi;	
								$info=~s/MULTISPECIES://gi;	
								$info=~s/multi-domain//gi;							
								$info=~s/domain//gi;
								
										
								$info=~s/subtype//gi;
								$info=~s/type//gi;
								$info=~s/Cas_//gi;
								$info=~s/Cas\///gi;
									
								$info=~s/^\s+//;
								$info=~s/\s+$//;		
								$info=~s/^,//g;
								$info=~s/,$//g;
								$info=~s/\/\s+/\//g;
								$info=~s/\s+/ /g;
								
								
								#print "$info\t$rec_line\n";	
								
								
								#--- now replace the info with proposed name and add the other info next to the proposed name
								#--- check if multiple proposed name is found suggesting fused Cas genes
								my %hash_of_proposed_names;
								my $proposed_name_found=0;
								foreach my $s_name(sort{length($b)<=>length($a)}keys %hash_of_synonyms_and_proposed_cas_name)
									{
										if($info=~/\b($s_name)\b/i)
											{
												my $m_name=$1;
												my $proposed_cas_name=$hash_of_synonyms_and_proposed_cas_name{$s_name};												
												$hash_of_proposed_names{$proposed_cas_name}=1;
												
												#---- replace the cas name with proposed cas_name
												$info=~s/$m_name/$proposed_cas_name/g;
												
												$proposed_name_found++;
											}
									}
								my $nof_proposed_cas_genes_found=keys %hash_of_proposed_names;
								
								#------------------------------- final cleaning up ---------
								
								$info=~s/-containing//;
								$info=~s/homolog//;
								$info=~s/^\s+//;
								$info=~s/\s+$//;
								
								$info=~s/\s+,\s+/,/;
								$info=~s/\s+\/\s+/\//;
								$info=~s/\s+\//\//;
								$info=~s/\/\s+/\//;
								
								$info=~s/-$//;
								$info=~s/^-//;
								
								$info=~s/\/$//;
								$info=~s/^\///;
								
								if($info=~/(\S+)\/(\S+)\/(\S+)/)
									{
										my $w1=$1;
										my $w2=$2;
										my $w3=$3;
										if($w1 eq $w2 and $w2 eq $w3)
											{
												$info=~s/$w1\/$w2\/$w3/$w1/;
											}	
									}
								if($info=~/(\S+)\/(\S+)/)
									{
										my $w1=$1;
										my $w2=$2;
										if($w1 eq $w2)
											{
												$info=~s/$w1\/$w2/$w1/;
											}	
									}									
									
								$info=~s/,/ /g;	
								$info=~s/-/ /g;	
								$info=~s/\s+/ /g;
								$info=~s/ramp/RAMP/g;
								$info=~s/SAG0894/Cas9/g;
								$info=~s/Cas3\'/Cas3/gi;
								#-----------------------------------------------------------
								if($info eq "Cas" or $info eq "cas")
									{
										$info="";
									}
								if(defined $info and $info=~/\S+/)
									{
										if($proposed_name_found>0)
											{
												$cas_family=$info;
											}
										else{
												$cas_family=$info." Unclassified_Cas_protein";
											}		
									}
								#else{
								#		
								#	}	
								
																
								if($nof_proposed_cas_genes_found>1)
									{
										#print "Multiple proposed_cas_gene_found in $rec_line\n\t$cas_family\n\n";
									}	
								#------	
									
								my $skip_t1=1;
								if($skip_t1==0)
								{	
								my @arr_t2=split('CRISPR-associated ',$rec_line);
								my $cas_info=$arr_t2[1]; $cas_info=~s/\"//g;
								my @arr_t3=split('\s+',$cas_info);	
								if(defined $arr_t3[1] and defined $arr_t3[1]=~/\S+/)
									{
										$arr_t3[1]=~s/,//g;
										$arr_t3[1]=~s/;//g;
										
										$arr_t3[1]=~s/-/_/g;
										
										if($arr_t3[0]=~/protein/i)
											{
												if($arr_t3[1]=~/\S+/)
													{
														$cas_family=$arr_t3[1]."-".$arr_t3[0];;
													}
												else{
														$cas_family="Unclassified_Cas_protein";
													}	
											}
										elsif($arr_t3[1]=~/family/i)
											{
												$cas_family=$arr_t3[0];
											}							
										else{
												$cas_family=$arr_t3[1]."-".$arr_t3[0];
											}	
									}	
								}	
								#----- get the cas_protein_accession
								
								my $j=$i-1;
								while($j>0)
									{
										if(defined $arr_info_line[$j] and $arr_info_line[$j]=~/protein_id\=/)
											{
												$cas_protein_acc=$arr_info_line[$j];chomp $cas_protein_acc; $cas_protein_acc=~s/\r//g; 
												$cas_protein_acc=~s/\/protein_id\=//; $cas_protein_acc=~s/^\s+//; $cas_protein_acc=~s/\"//g; $cas_protein_acc=~s/\.\d+$//;
												last;
											}
										if($arr_info_line[$j]=~/^     CDS/){last;}	
										$j--;
									}
									
																	
								#if($cas_protein_acc eq "NA")
								#	{
								#		print "\n\nError processing: $input_gbk_file\n\t$current_accession\tCRISPR\t$cas_gene_start\t$cas_gene_stop\t$cas_family\t$cas_protein_acc\n";exit;
								#	}	
								print TAB "$current_accession\tCRISPR\t$cas_gene_start\t$cas_gene_stop\t$cas_family\t$cas_protein_acc\n";
						}			
						
					
						
			}		
		close(TAB);
		
		
		#----- get the taxonomy (short path)
		
		my @arr_lines=`grep -E '^LOCUS|  ORGANISM' $input_gbk_file -A 10 | awk '{var i;if(\$0~/^LOCUS/){i=0;};if(\$0~/^DEFINITION/){i=1;}; if(\$0~/^  ORGANISM/){i=0;};if(\$0~/^REFERENCE/){i=1;}; if(i!=1){print \$0}}' >&1`;
		
		
		open(TAB,">>$tmp_dir\/$output_gene_position_file") or print "";
		flock(TAB,2);
		my $c_accession;
		for(my $i=0;$i<=$#arr_lines;$i++)
			{
				#chomp $line; $line=~s/\r//;
				
				
				#----- get the accession and strain [if present]
				if($arr_lines[$i]=~/^LOCUS    /)
					{
						my @arr_t1=split('\s+',$arr_lines[$i]);
						$c_accession=$arr_t1[1];
					}
				if($arr_lines[$i]=~/^  ORGANISM/)
					{						
						my $org=$arr_lines[$i]; $org=~s/^  ORGANISM//g; chomp $org; $org=~s/\r//g;$org=~s/^\s+//;
						my $tax_path="";
						my $j=$i+1;
						while($j<=$#arr_lines)
							{
								chomp $arr_lines[$j]; $arr_lines[$j]=~s/\r//g;$arr_lines[$j]=~s/^\s+//;
								$tax_path=$tax_path.$arr_lines[$j];
								
								if(defined $arr_lines[$j+1] and $arr_lines[$j+1]=~/^LOCUS    /)
									{
										$i=$j;
										last;
									}
								if($j==$#arr_lines)
									{
										$i=$j;
										last;
									}	
								$j++;	
							}
						$tax_path=$tax_path."; ".$org;	
						my @arr_t1=split(';',$tax_path);	
						my $species=$arr_t1[0];
						print TAB "$c_accession\tTAXONOMY\t0\t0\t$species\t$tax_path\n";
					}	
					
			}
			
		close(TAB);				
		
		
		return 1;
	}

sub process_fasta_file()
	{
		my($input_fasta_file,$arr_sequence_files,$hash_id_lookup_table)=@_;
		
		
		
		#----------- now open the user fasta file and extract the individual sequences, and pass them to first pilercr, then to CRISPRDetect -----
		open(RD,"$input_fasta_file") or print "$!: $input_fasta_file not found\n";


		my $seq_index=0;
		my $last_seq_id;
		my $seq="";
		while( my $line=<RD> )
			{				
				chomp $line; $line=~ s/\r//; $line=~ s/^\s+//; 
				#print "$line\n";
				
				if($line=~/^>/)
					{
						#print "matched\n";
						if($seq_index>0)
							{
								$seq=uc($seq);
								open(WR,">$tmp_dir\/$last_seq_id\.fna") or print "$!";
								flock(WR,2);
								print WR ">$last_seq_id\n";
								print WR "$seq\n";
								close(WR);
								
								push(@$arr_sequence_files,$last_seq_id);
								
								$last_seq_id="";
								$seq="";
							}						
						
						my $o_id=$line;	chomp $o_id;$o_id=~ s/\r//; 
						$o_id=~ s/^>//;
						#---- check if NCBI accession style ---
						if($o_id=~/^gi\|/)
							{
								my @arr_t1=split('\|',$o_id);  #>gi|851213623|ref|NZ_CP009517.1| Methanosarcina barkeri 3, complete genome
								$last_seq_id=$arr_t1[3];
								$last_seq_id=~s/\.\d+$//;								
								
							}
						elsif($o_id=~/^SEQUENCE_/)    #---- already processed and unique
							{
								$last_seq_id=$o_id;
							}	
						else{	
								#$o_id=~ s/[^a-zA-Z0-9_.-]/_/;
															
															
								my $time=&get_unique_id();
							
								my $t_o_id=$o_id; $t_o_id=~s/\s//g; $t_o_id=~s/[-_,;]//g; $t_o_id=~s/\W//g;  $t_o_id=uc($t_o_id);
								if(length($t_o_id)>50){$t_o_id=substr($t_o_id,-50);}
								$last_seq_id="SEQUENCE_".$seq_index.$time.$t_o_id;	
							}
						
						$hash_id_lookup_table->{$last_seq_id}=$o_id;
						
						#print WR ">$new_id\n";
												
						$seq_index++;
					}
				else{
						#---------- get the sequence ---
						#print "matched\n";
						
						$line=uc($line);
						$seq=$seq.$line;
					}																	 								 						 
				#print WR "$line\n";
												
			 }
		close(RD); 
		#------- now write the last sequence
		$seq=uc($seq);
		
		#print "\$last_seq_id=$last_seq_id\n";
		open(WR,">$tmp_dir\/$last_seq_id\.fna") or print "$!";
		flock(WR,2);
		print WR ">$last_seq_id\n";
		print WR "$seq\n";
		close(WR);
			 
		push(@$arr_sequence_files,$last_seq_id);
			 
			
		return 1;
	}
	


sub simplify_crt_formatted_arrays_and_append()
	{
		my($crt_output_file,$pilercr_output_file)=@_;
		
		#print "$crt_output_file,$pilercr_output_file,$out_put_dir\n";
		my $highest_array_index=0;
		my %pilercr_array_ranges;	
		
		my $skip=0;
		
		if($skip==1)
		{	
		#----- open the PILER-CR output file and get the array index --------
		open(RD,"$tmp_dir\/$pilercr_output_file") or print "$!";
		flock(RD,2);
		my @arr_pilercr=<RD>;
		close(RD);
		
		my $highest_array_index=0;
		my $acc;
		for(my $i=0;$i<$#arr_pilercr;$i++)
			{
				my $line=$arr_pilercr[$i];
				if($line=~/^Array/ and $arr_pilercr[$i+1]=~/^>/)
					{
						if($line=~/^Array (\d+)/)
							{
								if($1>$highest_array_index){$highest_array_index=$1;}
							}
						if($arr_pilercr[$i+1]=~/^>(\S+)/)
							{
								$acc=$1;
								chomp $acc; $acc=~s/\r+//g;
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
											$pilercr_array_ranges{$acc}{$range}=1;
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
		}
		
		
		#---- now read the arrays predicted by crt and convert them to piler-CR format --------------
		open(CRT,"$tmp_dir\/$crt_output_file") or print "\nThe File: $tmp_dir\/$crt_output_file not found \n";
		flock(CRT,2);
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
						#print "\$accession=$accession\t$line\t\$highest_array_index=$highest_array_index\n";#next;
						my $range;my $range_center; my $crt_array_length;
						
						if($arr_crt[$i]=~/Range: (\d+) - (\d+)/)
							{
								$range=$1."-".$2;
								$range_center=$1+int(($2-$1)/2);
								$crt_array_length=$2-$1+1;
							}	
						#--------- check if the array center belong to any piler-cr array already ----
						my $array_exist=0;
						if(defined $range)
							{
								foreach my $pilercr_range(keys %{$pilercr_array_ranges{$accession}})
									{
										#print "$pilercr_range\t$range_center<br>";
										my($r_start,$r_stop)=split('-',$pilercr_range);
										my $pilercr_array_length=$r_stop-$r_start+1;
										if($crt_array_length<$pilercr_array_length){$array_exist=1;last;}
										#if($range_center>=$r_start and $range_center<=$r_stop){$array_exist=1;last;}
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
								#print "\$arr_crt[$j]=$arr_crt[$j]";
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
										$heighest_repeat_occurrence_score=$hash_of_repeats{$repeat};
										$model_repeat=$repeat;
									}
							}
							
						#-------------- skip this array if no_of_repeats are less than 3-------------	
						#if($#arr_repeats<2){next;}
						
						#---- now skip the arrays below a certain degenaracy cutoff -----------------
						
							#--- check if there are at least 2 repeats with no degeneracy -----------
	
						#if($heighest_repeat_occurrence_score<2){next;}

						#--------- check if there is just one repeat in the array
						
						if($#arr_of_model_repeats>=0)
							{
								#----- now select the best model repeat by scoring the array
								my %tmp_hash1;
								foreach my $repeat(@arr_of_model_repeats)
									{

										my $array_degeneracy_score=&get_array_degeneracy_score($model_repeat,$repeat,\@arr_repeats);
										
										#print "\$array_degeneracy_score=$array_degeneracy_score\n";
										
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
						
						
						#----
						$model_repeat=~s/\s+$//;
						#print "\$model_repeat=$model_repeat\n";
					
						
						#----- now that the model repeat is found, convert the whole array on the go
						$highest_array_index++;
						open(WR,">>$tmp_dir\/$pilercr_output_file") or print "$!";
						flock(WR,2);
						
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
						print WR "\n\n";
						#print WR "SUMMARY BY SIMILARITY\n";
						close(WR);
						#--- end of current array ----
						$i=$k;
						#exit;
					}
			}
		
		
		return 1;
	}


sub get_current_array()
	{
		my($original_array,$current_array)=@_;
				
		#---------------------------------
			
		my $avg_spacer_length=0;		
		my $no_of_spacers=0;			
		my $index=0;	
		my $skip_line=0;	
		my $total_spacer_length=0;				
		for(my $k=0;$k<4;$k++)
			{
				chomp $$original_array[$k]; $$original_array[$k]=~ s/\r+//g;
				if($k<(2)){push(@{$current_array},$$original_array[$k]);}				
				elsif($k==(2)){push(@{$current_array},"Position\tRepeat\tSpacer\tComment");}
				elsif($k==(3)){push(@{$current_array},"========\t======\t======\t=======");}
			}	
		
		for(my $k=4;$k<=$#{$original_array};$k++)
			{
				#print "$arr_compiled_pilercr_output[$k]";
				chomp $$original_array[$k]; $$original_array[$k]=~ s/\r+//g;
				#push(@original_array,$array_cpr_lines[$k]);
				
				
				#print "$array_cpr_lines[$k]\n";
				#-------------------------------------------------------------------

				#if($k<(4)){next;}
				
				#---- now make a simpler array with 4 columns ------------------------------------
				my $current_line=$$original_array[$k]; 
				chomp $current_line; $current_line=~s/\r+//g;
				$current_line=~s/^\s+//;$current_line=~s/\s+/\t/g;
				
				if($current_line=~/^=/){push(@{$current_array},"========\t======\t======\t=======");$skip_line=1;}
				
				
				
				
				if($skip_line==0)   #---if skip_line is not set to skip, continue
				{	#next;}	
				
						
				#print "$current_line\n";
							
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq);	
				my $insertion_base_and_position="";
				
				if($tmp_array[0]==1 and $#tmp_array!=6)    #--- when there are no left flank (array starts at first base), check NC_014374: remember this is already fixed when I convert a CRT array to PILER-CR
					{
						$r_start=$tmp_array[0];
						$r_length=$tmp_array[1];
						$s_length=$tmp_array[3];
						$r_seq=$tmp_array[4];
						$s_seq=$tmp_array[5];
						
						$total_spacer_length=$total_spacer_length+length($s_seq);
						$no_of_spacers++;
					}		
				elsif($#tmp_array==6)
					{
						$r_start=$tmp_array[0];
						$r_length=$tmp_array[1];
						$s_length=$tmp_array[2];
						$r_seq=$tmp_array[5];
						$s_seq=$tmp_array[6];
						
						$total_spacer_length=$total_spacer_length+length($s_seq);
						$no_of_spacers++;
					}
				
				elsif($#tmp_array==5)
					{
						$r_start=$tmp_array[0];
						$r_length=$tmp_array[1];
						$s_length=0;
						$r_seq=$tmp_array[4];
						$s_seq="|";
					}
				elsif($#tmp_array==4)
					{
						$r_start=$tmp_array[0];
						$r_length=$tmp_array[1];
						$s_length=0;
						$r_seq=$tmp_array[4];
						$s_seq="|";
					}		
				#print "\n\n\$r_start=$r_start\n";		
				
				#print "\$r_start=$r_start  \$r_seq=$r_seq  \$s_seq=$s_seq\n";
				
				#--------- check if the repeats has gaps followed by a single dot (dot means a match) ---------
				if($r_seq=~/\-+\.$/) #----- right side end
					{
						$r_seq=~s/\.$/-/;
						my @arr_mr_t1=split('',$model_repeat);
						$s_seq=$arr_mr_t1[$#arr_mr_t1].$s_seq;
					}
				if($r_seq=~/^\.\-+/ and $k>4) #----- left side end
					{
						$r_seq=~s/^\./-/;
						my @arr_mr_t1=split('',$model_repeat);
						$s_seq=$arr_mr_t1[$#arr_mr_t1].$s_seq;
						
						#---- now fetch the previous record, and add the base to the spacer right
						my $previous_rec=$$current_array[$#{$current_array}];
						my ($pre_r_start,$pre_r_seq,$pre_s_seq,$tmp_1)=split('\t',$previous_rec);
						
						$pre_s_seq=$pre_s_seq.$arr_mr_t1[0];
						$$current_array[$#{$current_array}]="$pre_r_start\t$pre_r_seq\t$pre_s_seq\t";
					}
				
				#--------- check if the repeats has gaps followed by a single dot (dot means a match) ---------
				if($r_seq=~/(\w+)$/) #----- right side end
				{
						#print "\nR:$r_seq $s_seq\n";
						
						my $matching_word=$1;
						
						
						if(length($matching_word)>1 and $k<$#original_array-1)
							{
								my $gap_string="";
								for(my $i=0;$i<length($matching_word);$i++)
									{
										$gap_string=$gap_string."-";
									}
								
								$s_seq=$matching_word.$s_seq;
								
								$r_seq=~s/$matching_word$//;
								$r_seq=$r_seq.$gap_string;
								
								#print "2:$r_seq $s_seq\n";	
							}
						

					}

				if($r_seq=~/^(\w+)/) #----- right side end
					{
						#print "\nL:$r_seq $s_seq\n";
						my $matching_word=$1;
						
						if(length($matching_word)>1 and $k>4)
							{
								my $gap_string="";
								for(my $i=0;$i<length($matching_word);$i++)
									{
										$gap_string=$gap_string."-";
									}
								
								my $previous_rec=$$current_array[$#{$current_array}];
								my ($pre_r_start,$pre_r_seq,$pre_s_seq,$tmp_1)=split('\t',$previous_rec);
								
								$pre_s_seq=$pre_s_seq.$matching_word;
								$$current_array[$#{$current_array}]="$pre_r_start\t$pre_r_seq\t$pre_s_seq\t";
								
								$r_seq=~s/^$matching_word//;
								$r_seq=$gap_string.$r_seq;
								
								#print "2:$r_seq $s_seq\n";	
							}
						
					}
				
				
				
				
				#--------- check for no_of_gaps in the repeat ---------
				my $no_of_gaps=0;
					
				if($r_seq=~/-/)
					{
						$no_of_gaps= $r_seq=~ s/-/-/g;					
					}
				#---- address the coord difference
					
				if($#{$current_array}>3)
					{
						#print "\$#current_array=$#current_array\n";	
						my $last_record=$$current_array[$#{$current_array}];		
						#print "\n\$last_record=$last_record\n";	
													
						my($last_rec_start,$last_repeat,$last_spacer,$comment)=split('\t',$last_record);
								
						$last_repeat=~s/-//g;
						$last_spacer=~s/-//g;
								
						$r_start=$last_rec_start+length($last_repeat)+length($last_spacer);
						
						
						#$coord_diff=$coord_diff+$no_of_gaps;
						#if($first_occurrence_of_gap>0)
						#	{
						#		$r_start=$r_start-$coord_diff;
						#	}
						#else{$first_occurrence_of_gap++}		
					}
				
				#----- create a rec line and save it into current array --
				
				my $rec_line="$r_start\t$r_seq\t$s_seq\t";
				#print "$rec_line\n";
				push(@{$current_array},$rec_line);
				$index++;
				}
				
			} 
		
		
		if($no_of_spacers==0){print "@{$original_array}\n\n"; return(0,0);}#$pm->finish; exit;}
		else{	
			  $avg_spacer_length=int($total_spacer_length/$no_of_spacers);	
		    }
				
				
				
				
		#-------------------------------------------------
				
		return($avg_spacer_length,$total_spacer_length);
	}
			



sub check_array_for_misaligned_bases()
	{
		my($range,$accession,$model_repeat,$avg_spacer_length,$current_array,$modified_array)=@_;
		
		
		#---------- now open the sequence file and get the sequence string -------------
		open(SEQ,"$tmp_dir\/$accession\.fna") or print "Error 1: can't find $accession.fna $! <br>\n";
		my @arr_seq=<SEQ>;
		close(SEQ);
		my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;	
		#------------------------------------------------------------------------------------
		
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
				
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;							
				#print "\$current_line=$current_line\n";
							
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
						
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1];
				
				if(defined $tmp_array[2])
					{
						$s_seq=$tmp_array[2];
					}
				else{
						$s_seq="-";
					}		
				my $s_seq_1=$s_seq;$s_seq_1=~s/-//g; #if(not $s_seq_1 or $s_seq_1 eq ""){$s_seq_1="-";}
				
				
				#$s_seq=$right_flank;
				my $existing_insertion_bases="";
				my $existing_insertion_positions="";
				#my %hash_of_insertion_bases;
				#my $number_of_insertions=0;
				if(defined $tmp_array[3])
					{
						$comment=$tmp_array[3]; 
						
						if($tmp_array[3]!~/^\s{0,1}Del/)
							{
								#$number_of_insertions=&get_number_of_insertions($tmp_array[3],\%hash_of_insertion_bases);
								
								($existing_insertion_bases,$existing_insertion_positions)=split(' ',$comment);
								if($existing_insertion_positions)
									{
										$existing_insertion_positions=~s/\[//g;
										$existing_insertion_positions=~s/\]//g;
									}
							}	
					}
				else{
						$comment="";
					}	
																				
				$r_length=length($r_seq);
				#$s_length=length($s_seq);
				my $r_seq_t=&change_dots_to_bases($r_seq,$model_repeat);
				
				#----- logic: look for the longest dotted region, and split the $r_seq with that,
				#-----		   the longest dotted region can start from the first base, middle, or last base
				#-----			if the mismatch is on the right, add the whole repeat to the spacer and pass to clustalw to get the alignment with $model_repeat				
				
				
				#--------- check if the repeats has gaps followed by a single dot (dot means a match) ---------
				if($r_seq=~/\-+\.$/) #----- right side end
					{
						$r_seq=~s/\.$/-/;
						my @arr_mr_t1=split('',$model_repeat);
						$s_seq=$arr_mr_t1[$#arr_mr_t1].$s_seq;
					}
				if($r_seq=~/^\.\-+/ and $k1>4) #----- left side end
					{
						$r_seq=~s/^\./-/;
						my @arr_mr_t1=split('',$model_repeat);
						$s_seq=$arr_mr_t1[$#arr_mr_t1].$s_seq;
						
						#---- now fetch the previous record, and add the base to the spacer right
						my $previous_rec=$$current_array[$k1-1];
						my ($pre_r_start,$pre_r_seq,$pre_s_seq,$tmp_1)=split('\t',$previous_rec);
						
						$pre_s_seq=$pre_s_seq.$arr_mr_t1[0];
						$$current_array[$k1-1]="$pre_r_start\t$pre_r_seq\t$pre_s_seq\t";
					}
				
				my $tmp_repeat=$r_seq;$tmp_repeat=~s/\.//g;
				my $number_of_mismatch=length($tmp_repeat);
				
				if($r_seq=~/^(\.{0,})([^\.]\S+)/ and $number_of_mismatch>3 and $k1>4 )
					{
						#print "$r_seq\n";
						
						my $old_repeats_similarity_score=&get_similarity_score($r_seq_t,$model_repeat);
						
						
						
						
						
						my $tmp_r_seq=substr($r_seq,-$-[0]);
						#print "\nR:$r_seq \t $-[1] $+[1] $2 $tmp_r_seq\n";
								
						my $previous_rec=$$current_array[$k1-1];
						my ($pre_r_start,$pre_r_seq,$pre_s_seq,$tmp_1)=split('\t',$previous_rec);
						my $gapless_pre_s_seq=$pre_s_seq;$gapless_pre_s_seq=~s/-//g;						
						
						my $p_s_seq_start=$r_start-length($gapless_pre_s_seq);
						
						my $gapless_r_seq=$r_seq;$gapless_r_seq=~s/-//g;
						my $number_of_insertions=0;
								
						if(defined $comment and $comment ne "" and $comment!~/^Del/)
							{
								my %tmp_hash_of_insertions;
								$number_of_insertions=&get_number_of_insertions($comment,\%tmp_hash_of_insertions);
							}
						if(not defined $number_of_insertions){$number_of_insertions=0;}
						
						my $tmp_seq="";
						
						if(($p_s_seq_start-1)>=0 and ($p_s_seq_start-1)<length($species_seq))
							{
								$tmp_seq=substr($species_seq,$p_s_seq_start-1,length($gapless_pre_s_seq)+length($gapless_r_seq) + $number_of_insertions + length($s_seq_1)); # This one worked for the module fix_array_with_falsely_iden tified_repeats perfectly
							}
						
						#my $tmp_seq=substr($species_seq,$p_s_seq_start,length($gapless_pre_s_seq)+length($gapless_r_seq) + $number_of_insertions + length($s_seq_1)-1); # -1 for zero based coord
						
						
					
						my %hash_of_lines;		
						&run_clustalw_on_two_seqiuences($range,$accession,$tmp_seq,$model_repeat,\%hash_of_lines);
						
						#foreach my $index(sort{$a<=>$b}keys %hash_of_lines)
						#	{
						#		print "$hash_of_lines{$index}\n";
						#	}
							
						#----- get the previous spacer --------------------------
						my $spacer_repeat_spacer_line="";
						my $mr_line="";
						
						if(defined $hash_of_lines{0})
							{
								$spacer_repeat_spacer_line=$hash_of_lines{0};
							}
						if(defined $hash_of_lines{1})
							{		
								$mr_line=$hash_of_lines{1};
							}	
						
						my $no_of_leading_gaps=0;
						my $no_of_trailing_gaps=0;
						
						my $new_previous_spacer="";
						my $new_repeat_start="";
						my $new_repreat="";
						my $new_current_spacer="";
						my $new_comment="";
						
						if($mr_line=~/^(\-+)/)
							{
								$no_of_leading_gaps=length($1);
								
								$new_previous_spacer=substr($spacer_repeat_spacer_line,0,$no_of_leading_gaps);
								#$mr_line=~s/^$pre_s_seq//;
							}	
						if($mr_line=~/(\-+)$/)
							{
								$no_of_trailing_gaps=length($1);
								$new_current_spacer=substr($spacer_repeat_spacer_line,-$no_of_trailing_gaps);
							}		
						
						$new_repeat_start=$p_s_seq_start+length($new_previous_spacer);
								
								#--change $r_start accordingly ---
								#my $r_start_t1=$p_s_seq_start+length($pre_s_seq);
								
						my $aligned_repeat_region=substr($spacer_repeat_spacer_line,$no_of_leading_gaps,length($spacer_repeat_spacer_line)-$no_of_leading_gaps-$no_of_trailing_gaps);
						my $aligned_model_repeat_region=substr($mr_line,$no_of_leading_gaps,length($spacer_repeat_spacer_line)-$no_of_leading_gaps-$no_of_trailing_gaps);		
						
								
						#print "\$pre_s_seq=$pre_s_seq\t\$r_start=$r_start\t$aligned_repeat_region\t$new_current_spacer\n";
						
						
						#------------ now handle the inserted bases --------------------------------
						
						my @arr_top_seq=split('',$aligned_repeat_region);
						my @arr_bottom_seq=split('',$aligned_model_repeat_region);
						my $insertion_bases="";
						my $insertion_bases_positions="";
						
						for(my $i=0;$i<=$#arr_bottom_seq;$i++)
							{
								if($arr_bottom_seq[$i] eq "-")
									{
										if(defined $arr_bottom_seq[$i-1] and $arr_bottom_seq[$i-1] eq "-")
											{
												$insertion_bases=$insertion_bases.$arr_top_seq[$i];
												#$insertion_bases_positions=$insertion_bases_positions.$i;
											}
										else{	
												$insertion_bases=$insertion_bases.",".$arr_top_seq[$i];
												$insertion_bases_positions=$insertion_bases_positions.",".int($p_s_seq_start+$no_of_leading_gaps+$i);
											}
										$arr_top_seq[$i]="";
									}
							}
						$insertion_bases=~s/^,//;	
						$insertion_bases_positions=~s/^,//;
						
						$new_repreat=join("",@arr_top_seq);
						
						if($insertion_bases ne "")
							{
								$new_comment="$insertion_bases [$insertion_bases_positions]";	
							}	
						
						#--- now check if this alignment is better than previous one
						my $new_repeats_similarity_score=&get_similarity_score($new_repreat,$model_repeat);
						
						#print "$new_repeats_similarity_score\t$old_repeats_similarity_score\t\$r_seq=$r_seq\t\$new_comment=$new_comment\n";
						
						if($new_repeats_similarity_score>$old_repeats_similarity_score)
							{
								
								#--- first change the previous record ---
								$$current_array[$k1-1]="$pre_r_start\t$pre_r_seq\t$new_previous_spacer\t$tmp_1";
								
								#--- now change the current record ------
								$new_repreat=&change_bases_to_dots($new_repreat,$model_repeat);
								
								$r_start=$new_repeat_start;
								$r_seq=$new_repreat;
								$s_seq=$new_current_spacer;
								$comment=$new_comment;
							}
									
					}
				
				
				#--------- check if the repeats has mismatches at the right side end ---------
				my $skip1=0;
				if($skip1==1)
				{
				
				if($r_seq=~/(\w+)\.{0,}$/) #----- right side end
				{
						my $matching_word=$1;
						
						
						
						
						#--- check the right side part where the first mismatch occured
						if($r_seq=~/^(\.{0,})([^\.]\S+)/)
							{
								my $tmp_r_seq=substr($r_seq,-$-[0]);
								#print "\nR:$r_seq \t $-[1] $+[1] $2 $tmp_r_seq\n";
								
								
								
							}
						
						
						
						
						
						if(length($matching_word)>1 and $k1< ($#{$current_array}-1))
							{
								my $gap_string="";
								for(my $i=0;$i<length($matching_word);$i++)
									{
										$gap_string=$gap_string."-";
									}
								
								$s_seq=$matching_word.$s_seq;
								
								$r_seq=~s/$matching_word$//;
								$r_seq=$r_seq.$gap_string;
								
								#print "2:$r_seq $s_seq\n";	
							}
						

					}

				if($r_seq=~/^(\w+)/) #----- left side end
					{
						#print "\nL:$r_seq $s_seq\n";
						my $matching_word=$1;
						
						if(length($matching_word)>1 and $k1>4)
							{
								my $gap_string="";
								for(my $i=0;$i<length($matching_word);$i++)
									{
										$gap_string=$gap_string."-";
									}
								
								my $previous_rec=$current_array[$k1-1];
								my ($pre_r_start,$pre_r_seq,$pre_s_seq,$tmp_1)=split('\t',$previous_rec);
								
								$pre_s_seq=$pre_s_seq.$matching_word;
								$current_array[$k1-1]="$pre_r_start\t$pre_r_seq\t$pre_s_seq\t";
								
								$r_seq=~s/^$matching_word//;
								$r_seq=$gap_string.$r_seq;
								
								#print "2:$r_seq $s_seq\n";	
							}
						
					}
				}
				
				$s_seq=~s/-+$/-/;
				$$current_array[$k1]="$r_start\t$r_seq\t$s_seq\t$comment";
			}
			
		return 1;
	}



sub fix_array_with_falsely_identified_repeats()
	{
		my($range,$accession,$model_repeat,$avg_spacer_length,$current_array,$modified_array)=@_;
		
		
		#---------- now open the sequence file and get the sequence string -------------
		open(SEQ,"$tmp_dir\/$accession\.fna") or print "Error 1: can't find $accession.fna $! <br>\n";
		my @arr_seq=<SEQ>;
		close(SEQ);
		my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;	
		#------------------------------------------------------------------------------------
		
		for(my $m=3;$m>=0;$m--)
			{
				unshift(@{$modified_array},$$current_array[$m]);
								#shift(@modified_array);
			}
			
		
		
		
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++) #---- starts from 5
			{
				
				
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;							
				#print "C=$k1: $current_line\n";
							
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
						
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1];
				
				if(defined $tmp_array[2])
					{
						$s_seq=$tmp_array[2];
					}
				else{
						$s_seq="-";
					}		
				my $s_seq_1=$s_seq;$s_seq_1=~s/-//g; #if(not $s_seq_1 or $s_seq_1 eq ""){$s_seq_1="-";}
				
				
				#$s_seq=$right_flank;
				my $existing_insertion_bases="";
				my $existing_insertion_positions="";
				#my %hash_of_insertion_bases;
				#my $number_of_insertions=0;
				if(defined $tmp_array[3])
					{
						$comment=$tmp_array[3]; 						
					}
				else{
						$comment="";
					}	
																				
				$r_length=length($r_seq);
				#$s_length=length($s_seq);
				my $r_seq_t=&change_dots_to_bases($r_seq,$model_repeat);
				
				#----- logic: look for the longest dotted region, and split the $r_seq with that,
				#-----		   the longest dotted region can start from the first base, middle, or last base
				#-----			if the mismatch is on the right, add the whole repeat to the spacer and pass to clustalw to get the alignment with $model_repeat				
				
				
				#--------- check if the repeats has gaps followed by a single dot (dot means a match) ---------

				
				my $tmp_repeat=$r_seq;$tmp_repeat=~s/\.//g;
				my $number_of_mismatch=length($tmp_repeat);
				
				
				if($r_seq=~/^(\.{0,})([^\.]\S+)/ and $number_of_mismatch>3 and $k1>4 )
					{
						#print "$r_seq\n";
						
						my $old_repeats_similarity_score=&get_similarity_score($r_seq_t,$model_repeat);
						
						
						
						
						
						my $tmp_r_seq=substr($r_seq,-$-[0]);
						#print "\nR:$r_seq \t $-[1] $+[1] $2 $tmp_r_seq\n";
								
						my $previous_rec=$$modified_array[$#{$modified_array}];
						my ($pre_r_start,$pre_r_seq,$pre_s_seq,$tmp_1)=split('\t',$previous_rec);
						my $gapless_pre_s_seq=$pre_s_seq;$gapless_pre_s_seq=~s/-//g;						
						
						my $p_s_seq_start=$r_start-length($gapless_pre_s_seq);
						
						my $gapless_r_seq=$r_seq;$gapless_r_seq=~s/-//g;
						my $number_of_insertions=0;
								
						if(defined $comment and $comment ne "" and $comment!~/^Del/)
							{
								my %tmp_hash_of_insertions;
								$number_of_insertions=&get_number_of_insertions($comment,\%tmp_hash_of_insertions);
							}
						if(not defined $number_of_insertions){$number_of_insertions=0;}
						
						my $tmp_seq="";
						
						if(($p_s_seq_start-1)>=0 and ($p_s_seq_start-1)<length($species_seq))
							{
								$tmp_seq=substr($species_seq,$p_s_seq_start-1,length($gapless_pre_s_seq)+length($gapless_r_seq) + $number_of_insertions + length($s_seq_1)); # -1 for zero based coord
							}	
						
						
						
						#---- if similarity score < 35%  and previous spacer length < median_spacer_length:  add the repeat to previous spacer and delete this row ----
						
						my $median_spacer_length=&get_median_spacer_length(\@{$current_array});
						
						if($old_repeats_similarity_score < int(length($model_repeat)/2) and length($gapless_pre_s_seq)<$median_spacer_length)
							{
								#print "\$old_repeats_similarity_score=$old_repeats_similarity_score and \$median_spacer_length=$median_spacer_length\n";
								#print "$current_line\n";
								
								#--- modify the previous spacer, and delete the current record --------------------	
								my $new_p_rec_line="$pre_r_start\t$pre_r_seq\t$tmp_seq\t$tmp_1";	
								#print "M=$#{$modified_array}: $new_p_rec_line\n\n";						
								$$modified_array[$#{$modified_array}]=$new_p_rec_line;
								
								#push(@{$modified_array},$$current_array[$#{$current_array}]);
								#delete $$current_array[$k1];
								#last;
								#next;
								
							}
						else{
								#print "P=$k1 $$current_array[$k1]\n\n";
								push(@{$modified_array},$$current_array[$k1]);
							}
						
						
						
						
									
					}
				else{		
						
						push(@{$modified_array},$$current_array[$k1]);
						#print "P=$#{$modified_array}: $$modified_array[$#{$modified_array}]\n";
					}
				

				#$s_seq=~s/-+$/-/;
				#$$current_array[$k1]="$r_start\t$r_seq\t$s_seq\t$comment";
			}
		
		push(@{$modified_array},$$current_array[$#{$current_array}]);	
		
		#----- now undef current array and assign modified_array to current_array
		#undef $current_array;
		@{$current_array}=@{$modified_array};
		
		#foreach my $index (0 .. $#{$current_array}) 
		#	{
		#		delete $$current_array[$index];
		#	}
		#foreach my $index1 (0 .. $#{$modified_array}) 
		#	{
		#		#print "$$modified_array[$index1]\n";
		#		$$current_array[$index1]=$$modified_array[$index1];
		#	}	
			
		return 1;
	}



sub get_number_of_insertions()
	{
		my($comment,$hash_of_insertion_positions)=@_;
		
		my $number_of_insertions=0;
		

		
		if($comment!~/^Del/)
			{
								
				my @tmp_arr1=split(' ',$comment);
														
														#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
				my $insertion_bases=$tmp_arr1[0];
				my $insertion_positions=$tmp_arr1[1];
								
				if($insertion_positions and $insertion_positions!~/\D/)
					{
								
					   $insertion_positions=~s/\[//g;
					   $insertion_positions=~s/\]//g;
																#print "\$cur_insertion_bases=$cur_insertion_bases\n";
												
						my @tmp_arr2=split(',',$insertion_bases);
						my @tmp_arr3=split(',',$insertion_positions);
									
						for( my $p=0;$p<=$#tmp_arr3;$p++)
							{
								my $pos=$tmp_arr3[$p];
								if(defined $tmp_arr2[$p])
									{
										chomp $tmp_arr2[$p];$tmp_arr2[$p]=~s/^\s+//;
										if($tmp_arr2[$p] ne "" and $tmp_arr2[$p] !~/\w+/)
											{
												$hash_of_insertion_positions->{$pos}=$tmp_arr2[$p];
											}
									}									
										#print "A:\$hash_of_insertion_positions{$p}=$hash_of_insertion_positions{$p}\n";
							}	
																		
						$insertion_bases=~s/,//g;
						$no_of_insertions=length($insertion_bases);						

					}

			}
		if(not defined $no_of_insertions or $no_of_insertions eq ""){$no_of_insertions=0;}
		return $no_of_insertions;				
	}

sub run_clustalw_on_two_seqiuences()
	{
		my($range,$accession,$seq1,$model_repeat,$hash_of_lines)=@_;
		
		$arr_spacers[0]=$seq1;
		$arr_spacers[1]=$model_repeat;
		#-----------------------------------------------------------------------------------------------------------------------------------------------------	
	
		my $time_1 = &get_unique_id();	   
		
		my $tmp_inputfile=	$accession.$range.$time_1."_tmp_spacers.txt";
		my $tmp_output=		$accession.$range.$time_1."_tmp_spacers.aln";
		my $tmp_output_dnd=	$accession.$range.$time_1."_tmp_spacers.dnd";
		
		
		
		
		open(FA,">$tmp_dir\/$tmp_inputfile") or print "$!";
		flock(FA,2);
		my $spacer_index=0;
		foreach my $spacer(@arr_spacers)
			{
				print FA ">S_$spacer_index\n$spacer\n";
				$spacer_index++;
			}
		close(FA);
		
		#---- run clustalW to get the repeat and spacers -----------------------------------------------------------------------------------------------------
		my @arr_output;
		#print "clustalw -INFILE=$tmp_dir\/$tmp_file -OUTFILE=$tmp_dir\/$tmp_output -QUIET -ALIGN -GAPOPEN=0.01 -GAPEXT=0.01 -OUTORDER=INPUT\n";
		system("clustalw -INFILE=$tmp_dir\/$tmp_inputfile -OUTFILE=$tmp_dir\/$tmp_output -QUIET -ALIGN -NUMITER=50 -NOWEIGHTS -ENDGAPS -OUTORDER=INPUT >/dev/null 2>&1");
		
		if(-e "$tmp_dir\/$tmp_output")
			{
				open(RD,"$tmp_dir\/$tmp_output");
				@arr_output=<RD>;
				close(RD);
				
				unlink("$tmp_dir\/$tmp_inputfile");	
				unlink("$tmp_dir\/$tmp_output");
				unlink("$tmp_dir\/$tmp_output_dnd");
			}
		#------------------------------------------------------------------------------
		for(my $i=0;$i<=$#arr_output;$i++)
			{
				my $line=$arr_output[$i];
				chomp $line; $line=~s/\r$//;
				
				if($line=~/^S_/)
					{	
											
						
						$line=~s/^S_//; $line=~s/\s+/\t/g;
						my @arr_t1=split('\t',$line);
						
						if(defined $hash_of_lines->{$arr_t1[0]})
							{
								$hash_of_lines->{$arr_t1[0]}=$hash_of_lines->{$arr_t1[0]}.$arr_t1[1];
							}
						else{								
								$hash_of_lines->{$arr_t1[0]}=$arr_t1[1];
							}
						#print "$line\n";		
					}
				elsif($line=~/^\s{16}/)
					{
						$line=~s/^\s{16}//;
						
						my $hash_index=$#arr_spacers+1;
						if($arr_output[$i-1]=~/^S_(\d+)/){$hash_index=$1+1;}
						
						if(defined $hash_of_lines->{$hash_index})
							{
								$hash_of_lines->{$hash_index}=$hash_of_lines->{$hash_index}.$line;
							}
						else{								
								if($line ne "")
									{
										$hash_of_lines->{$hash_index}=$line;
										#print "$line\n";
									}	
							}
					}	
					
			}

		
		return 1;	
	}




sub fix_arrays_with_insertion_in_repeat()
	{
		
		#print "\n";
		my($range,$accession,$model_repeat,$avg_spacer_length,$current_array,$modified_array)=@_;
		
		#print "Going to handle insertion(s)  for $accession :\n";
		
		
		my $case_found=0;
		
		if($model_repeat!~/-/){return($model_repeat,$case_found);}
		
		# cases:
		#1. When there are gap(s) at the end of model repeat, reduce the repeat length to avoid the gap(s), add the bases (if present) to the flank or spacers
		#2. 
		
		#		>gi|219883096|ref|NC_011880.1| Cyanothece sp. PCC 7425 plasmid pP742501, complete sequence
		#
		#       Pos  Repeat     %id  Spacer  Left flank    Repeat                                       Spacer
		#==========  ======  ======  ======  ==========    =========================================    ======
		#     98168      41    75.6      37  TTCAAACAGT    -----------..C.........................--    ATCCGAGACAGCGGCCGCGGCAGCGGGTTAGGCTCTC
		#     98246      41    97.6      32  TTAGGCTCTC    ---...............................A.....-    CGGAGGAGTCTTCGTCCCGATCACCGAACATG
		#     98319      41    90.2          ACCGAACATG    ACA.....................................A    CAGGTACCAG
		#==========  ======  ======  ======  ==========    =========================================
		#         3      41              34                ---GTCCCTACTTGTTAGGGAAACCAATTGAATGGAAACT-
		
		#print "LF:$left_flank\nRF:$right_flank\n";
		
		#------ leading gaps ----------------------------
		my $no_of_leading_gaps_in_mr=0;						
		if($model_repeat=~ /^(\-+)/)
			{
				#my $gap_string=$1;												
				#print "\$1=$1\n";
				$no_of_leading_gaps_in_mr =length($1);								
				#print "\$1=$1\t$no_of_leading_gaps_in_mr\n";
			}	
		
		
		#---- trailing gaps ------------------------------	
		my $no_of_trailing_gaps_in_mr=0;			
		if($model_repeat=~ /(\-+)$/)
			{
				#my $gap_string=$1;											
				#print "\$1=$1\n";
				$no_of_trailing_gaps_in_mr =length($1);
								
				#print "\$1=$1\t$no_of_trailing_gaps_in_mr\n";
			}
		
		#--- intermediate gaps ---------------------------
		my $no_of_gaps_in_mr=0;
						
		if($model_repeat=~ /(\-+)/)
			{
				#my $gap_string=$1;									
				#print "\$1=$1\n";
				$no_of_gaps_in_mr =length($1);
								
				#print "\$1=$1\t$no_of_leading_gaps_in_mr\n";
			}
		
		my $coord_diff=0;
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
				#print "@{$current_array-[0]}\n";
				my $current_line=$$current_array[$k1]; 
				chomp $current_line; $current_line=~s/\r+//g;
					$current_line=~s/^\s+//;#$current_line=~s/\s+/\s/g;
							
					#print "$current_line\n";
							
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq);
				my $comment="";
				
				my $insertion_base_and_position="";
				
						
						#if($#tmp_array==6)
						#	{
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1];
				
				if($tmp_array[2])
					{
						$s_seq=$tmp_array[2];
					}
				else{
						$s_seq="-";
					}	
					
				my %tmp_hash_of_insertions;	
				my $existing_insertion_bases="";
				my $existing_insertion_positions="";
				if(defined $tmp_array[3] )
					{
						
						$comment=$tmp_array[3];
						
						#print "\$comment=$comment\t\$tmp_array[3]=$tmp_array[3]\n";
						
						if(defined $comment and $comment ne "" and $comment!~/^Del/)
							{
								#my %tmp_hash_of_insertions;	
								$number_of_insertions=&get_number_of_insertions($comment,\%tmp_hash_of_insertions);
								
								#foreach my $pos(keys %tmp_hash_of_insertions)
								#	{
								#		print "$pos->$tmp_hash_of_insertions{$pos}\n";
								#	}
							}
						
						
						if($tmp_array[3]!~/^\s{0,1}Del/)
							{
								($existing_insertion_bases,$existing_insertion_positions)=split(' ',$comment);
								if($existing_insertion_positions)
									{
										$existing_insertion_positions=~s/\[//g;
										$existing_insertion_positions=~s/\]//g;
									}
							}
							
								
					}			
				$r_length=length($r_seq);
				$s_length=length($s_seq);
				
				#-----------------------------------------------------
				my($pre_r_start,$pre_r_seq,$pre_s_seq);
				my $pre_comment="";	
				if($k1>4) # get the previous record, as it will be used in many occations
					{
						my @arr_previous_rec=split('\t',$$current_array[$k1-1]);
						$pre_r_start=$arr_previous_rec[0];
						$pre_r_seq=$arr_previous_rec[1];
						
						if($arr_previous_rec[2])
							{
								$pre_s_seq=$arr_previous_rec[2];
							}
						else{
								$pre_s_seq="-";
							}	
								
						if($arr_previous_rec[3])
							{
								$pre_comment=$arr_previous_rec[3];
							}	
					}				
				
				#-------- now split the repeat in an array, which will be used through out the sub
				my @array_repeat=split('',$r_seq);	
				my @array_model_repeat=split('',$model_repeat);
				
				
				
				
				#-------- first check left side---------------------------------------------------------------------------------
				if($model_repeat=~/^-/)
					{
						#--- chop and add the leading base(s) from each repeat and add to the left flank/previous spacer sequence
						
						#-- case 1: if the insertion base is a gap, delete the gap
						
							
						for(my $i=0;$i<$no_of_leading_gaps_in_mr;$i++)
							{
								#if($array_repeat[$i] eq "-")
								#	{
								#		$array_repeat[$i]="";
								#		#shift(@array_repeat);
								#		$r_seq=join("",@array_repeat);
								#				
								#		#$r_start=$r_start-$coord_diff;
								#	}
								
								#---- check if there are any insertion bases at this coord
								my $current_pos=$r_start+$i;
								my $current_insertion_bases="";
								if(defined $tmp_hash_of_insertions{$current_pos})
									{
										$current_insertion_bases=$tmp_hash_of_insertions{$current_pos};
										#print "\$tmp_hash_of_insertions{$current_pos}=$tmp_hash_of_insertions{$current_pos}\n";
										delete $tmp_hash_of_insertions{$current_pos};
									}
									
								if($array_repeat[$i] ne "-")
									{
										#$coord_diff=$coord_diff+1;	
										
										if($k1==4)#--- check if this is the first record or not, if first record, add the base to the left flank
											{
												#$left_flank=$left_flank.$array_repeat[$i];												
												
												#$array_repeat[$i]="";												
												#$r_seq=join("",@array_repeat);
												
												$r_start=$r_start-1;																							
											}
										else{
												#---- add the base to the previous spacer sequence
												
												$pre_s_seq=$pre_s_seq.$current_insertion_bases.$array_repeat[$i];
												#$array_repeat[$i]="";$r_seq=join("",@array_repeat);
												
												$r_start=$r_start-1;
											}	
									}
								elsif($current_insertion_bases ne "")
									{
										$pre_s_seq=$pre_s_seq.$current_insertion_bases;
									}
									
								$array_repeat[$i]="";
								#shift(@array_repeat);
								$r_seq=join("",@array_repeat);	
									
							}	
					}
			
				
				#------------- now check right side -----------------------------------------------------------------------------
				if($model_repeat=~/-$/)
					{
						#--- chop and add the leading base(s) from each repeat and add to the left flank/previous spacer sequence
						
						#-- case 1: if the insertion base is a gap, delete the gap
						
							
						for(my $i=$#array_repeat;$i>($#array_repeat-$no_of_trailing_gaps_in_mr);$i--)
							{
								#---- check if there are any insertion bases at this coord
								my $current_pos=$r_start+$i;
								my $current_insertion_bases="";
								if(defined $tmp_hash_of_insertions{$current_pos})
									{
										$current_insertion_bases=$tmp_hash_of_insertions{$current_pos};
										#print "\$tmp_hash_of_insertions{$current_pos}=$tmp_hash_of_insertions{$current_pos}\n";
										delete $tmp_hash_of_insertions{$current_pos};
									}
								if($array_repeat[$i] ne "-")
									{
										#$coord_diff=$coord_diff+1;	
										
										if($k1==($#{$current_array}-1))#--- check if this is the last record or not, 
											{
												#$right_flank=$array_repeat[$i].$right_flank;												
												
												#$array_repeat[$i]="";												
												#$r_seq=join("",@array_repeat);
												
												#$r_start=$r_start-1;																							
											}
										else{
												#---- add the base to the previous spacer sequence
												
												$s_seq=$current_insertion_bases.$array_repeat[$i].$s_seq;
												#$array_repeat[$i]="";$r_seq=join("",@array_repeat);
												
												#$r_start=$r_start-1;
											}	
									}
								elsif($current_insertion_bases ne "")
									{
										$s_seq=$current_insertion_bases.$s_seq;
									}	
								
									
								$array_repeat[$i]="";
								#shift(@array_repeat);
								$r_seq=join("",@array_repeat);	
							}	
					}
			
				#------------- now check middle  -----------------------------------------------------------------------------
				
				if($model_repeat=~/\S-\S/)
					{
						my $insertion_string="";
						my $insertion_position="";
						
						#--- chop and add the leading base(s) from each repeat and add to the left flank/previous spacer sequence
						
						#-- case 1: if the insertion base is a gap, delete the gap
						
							
						for(my $i=$no_of_leading_gaps_in_mr;$i<=($#array_repeat-$no_of_trailing_gaps_in_mr);$i++)
							{
								if($array_model_repeat[$i] eq "-")  #------------------- check model repeat bases 
									{
										if(defined $array_model_repeat[$i-1] and $array_model_repeat[$i-1] eq "-")
											{
												if($array_repeat[$i] ne "-")
													{
														$insertion_string=$insertion_string.$array_repeat[$i];
													}
												#my $current_pos=$r_start+$i;
												#$insertion_position=$insertion_position.",".$current_pos;
											}
										else{
												if($array_repeat[$i] ne "-")
													{
														$insertion_string=$insertion_string.",".$array_repeat[$i];
														my $current_pos=$r_start+$i;
														
														$insertion_position=$insertion_position.",".$current_pos;
													}
												else{    #--- this is bug fix, where the first base of insertion is '-', example:  .-GTAGAAAC...................... in NC_007645's one of the arrays
														
														#$insertion_string=$insertion_string.",".$array_repeat[$i];
														my $current_pos=$r_start+$i;														
														$insertion_position=$insertion_position.",".$current_pos;
													}		
											}	
											
										$array_repeat[$i]="";
										#shift(@array_repeat);
										#$r_seq=join("",@array_repeat);			
									}							
							}
						$r_seq=join("",@array_repeat);	
						
						$insertion_string=~s/^,//;	$insertion_string=~s/^-+//g;chomp $insertion_string;
						$insertion_position=~s/^,//; chomp $insertion_position;
						
						#print "\$insertion_string=$insertion_string \t \$insertion_position=$insertion_position\n";
						#--------------------------------------------------------------------------------------
						if($insertion_string)  #--- do not check existance of $insertion_position for this purpose
							{
								
								############################### the following block is not applicable here, as this is the first time when program will create insertion bases
								
								#my $in_positions;
								#if($insertion_position=~/,/)
								#	{
								#		my @arr_t_1=split($insertion_position);
								#		my @arr_t_2;
								#		foreach my $p(@arr_t_1)
								#			{
								#				my $t_pos=$r_start+$p;
								#				push(@arr_t_2,$t_pos);
								#			}
								#		$insertion_position=join(",",@arr_t_2);	
								#	}
								#else{	
								#		#if($insertion_position!~/^[0-9]/)
								#		#{
								#		#	print "\$current_line=$current_line\n\$insertion_position=$insertion_position \$insertion_string=$insertion_string \$accession=$accession\n";
								#		#}									
								#		$insertion_position=$r_start+$insertion_position;
								#	}
								
								
														
								if($existing_insertion_bases ne "" and defined $existing_insertion_positions)
									{
										$insertion_string=$insertion_string.",".$existing_insertion_bases;
										
										$insertion_position=$insertion_position.",".$existing_insertion_positions;
									}
																						
								$comment= $insertion_string." [$insertion_position]";		
							}
						elsif($existing_insertion_bases ne "" and defined $existing_insertion_positions)
							{
								$comment="$existing_insertion_bases [$existing_insertion_positions]";
							}		
					}
			
			
				#----- now re-create the comment -------------------------------
				
				if(keys %tmp_hash_of_insertions > 0)
				{
				my @arr_t1;
				my @arr_t2;
				foreach my $pos(keys %tmp_hash_of_insertions)
					{
						push(@arr_t1,$pos);
						push(@arr_t2,$tmp_hash_of_insertions{$pos});
						#print "$pos->$tmp_hash_of_insertions{$pos}\n";
					}
					
					$comment=join(",",@arr_t2)."[".join(",",@arr_t1)."]";
				}
				else{$comment="";}
				#----- now update the current and old record
				
				#print "$pre_r_start\t$pre_r_seq\t$pre_s_seq\t$pre_comment\n";
				#print "$r_start\t$r_seq\t$s_seq\t$comment\n\n";
				if($k1>4)
					{
						if(not defined $pre_s_seq or $pre_s_seq eq ""){$pre_s_seq="-";}	
						$$current_array[$k1-1]="$pre_r_start\t$pre_r_seq\t$pre_s_seq\t$pre_comment";
					}
					
				if(not defined $s_seq or $s_seq eq ""){$s_seq="-";}		
				$$current_array[$k1]="$r_start\t$r_seq\t$s_seq\t$comment";
					
			}#end of for $k1
		#print "\nLF:$left_flank\nRF:$right_flank\n";
			
			
			
		#---- at the end remove all the gap(s) from the model_repeat seq		
		if($model_repeat=~/-/)
			{
				$model_repeat=~s/-//g;	
				$case_found=1;
			}
		return($model_repeat,$case_found);
	}





sub extend_array()
	{		
		
		my($range,$max_gap_between_crisprs,$ea_dynamic_search,$allowed_percent_similarity,$accession,$model_repeat,$avg_spacer_length,$current_array,$modified_array)=@_;
		
		#print "Going to search for extension in $accession:\n\n";
		
		my $case_found=0;
		my $left_flank="";
		my $right_flank="";
		#---------- now open the sequence file and get the sequence string -------------
		open(SEQ,"$tmp_dir\/$accession\.fna") or print "Error 1: can't find $accession.fna $! <br>\n";
		my @arr_seq=<SEQ>;
		close(SEQ);
		my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;		
		#---- logic ------------------------------------------------------------------------
		#- Step 1: 
		

		my $devision_length=int(length($model_repeat)/4);
		#my $k1=$#{$current_array}-1;
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
				
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;							
				#print "\$current_line=$current_line\n";
							
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
						
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1];
				$s_seq=$tmp_array[2];
				
				#$s_seq=$right_flank;
				my $existing_insertion_bases="";
				my $existing_insertion_positions="";
				if($tmp_array[3])
					{
						$comment=$tmp_array[3];
						
						if($tmp_array[3]!~/^\s{0,1}Del/)
							{
								($existing_insertion_bases,$existing_insertion_positions)=split(' ',$comment);
								if($existing_insertion_positions)
									{
										$existing_insertion_positions=~s/\[//g;
										$existing_insertion_positions=~s/\]//g;
									}
							}	
					}
				else{
						$comment="";
					}	
																				
				$r_length=length($r_seq);
				#$s_length=length($s_seq);
				
				
				#------------- left flank -------------------
				if($k1==4)
					{
						push(@{$modified_array},$current_line);
						#------------------------now search for degenerated repeats and spacer
						#next;
					
						#print "Inside\n\n\n";
						
						
						
						my @arr_previous_rec=split('\t',$$modified_array[0]);		#split('\t',$$current_array[$k1-1]);												
						
						my $last_repeat_start=$arr_previous_rec[0];
						my $last_repeat_seq=$arr_previous_rec[1];
						
						my $last_spacer_seq;
						
						if(defined $arr_previous_rec[2])
							{
								$last_spacer_seq=$arr_previous_rec[2];  
							}
						else{
								$last_spacer_seq="";
							}	
						
						#--- use while condition -----------------------------------------------------------
						my ($original_best_fitting_pattern,$best_fitting_pattern,$insertion_base,$position);
						my $pos;
						$best_fitting_pattern="-";
						
						
						my $skip_this_rec=0;
						my $extension_found=0;
						
						
						my $potential_previous_repeat="";
						my $potential_previous_spacer="";
						my $no_of_tries=0;
						while($best_fitting_pattern ne "")
							{
								
								$no_of_tries++;
								if($no_of_tries>5){last;}
								#print "$$modified_array[0]\n";
								
								
								
								#-------- find the repeat_start position from the last record -----------
								my @arr_cur_rec=split('\t',$$modified_array[0]);	
								
								#print "\t@arr_pre_rec\n";
								my $cur_r_start=$arr_cur_rec[0];
								my $cur_r_seq=$arr_cur_rec[1];	my $gapless_cur_r_seq=$cur_r_seq; $gapless_cur_r_seq=~s/-//g;										
								my $cur_s_seq="";
								if($arr_cur_rec[2])
									{
										$cur_s_seq=$arr_cur_rec[2]; $cur_s_seq=~s/-//g;
									}
								
								my $cur_comment=$arr_cur_rec[3];
								
								my $cur_no_of_insertions=0;
								my ($cur_insertion_bases,$cur_insertion_positions);
								#if($arr_pre_rec[3])
								#	{
								#		#$comment=$tmp_array[3];
								#		($pre_insertion_bases,$pre_insertion_positions)=split(' ',$pre_comment);										
								#		$pre_insertion_bases=~s/,//g;
								#		$pre_no_of_insertions=length($pre_insertion_bases);
								#	}	
								
														
								my $string_length=0;
								
								if(length($last_spacer_seq)>15 and length($last_spacer_seq)<70)
									{
										$string_length=int((length($last_repeat_seq) + length($last_spacer_seq)*1.33*$no_of_tries));					#	-length($potential_previous_spacer);								
									}
								else{
										$string_length=int(length($last_repeat_seq)*2.5);
									}		
								
								if($cur_r_start-$string_length-1 <0 )
									{
										
										$potential_previous_spacer="";
										#-------- make a dummy string uning remaining bases and '-'s
										
										if(($cur_r_start-1) > 0)
											{
												$potential_previous_spacer= substr($species_seq,0,$cur_r_start-1);
												#print "**** $r_start-1 : $previous_spacer_bases_to_chop\n\n";
											}	
										
										
										#-------- make a dummy string and take substr of it : currently just ending the process----------------- Example: NC_019776
										#print "\n\n############## here #################\n\n";
										my $dummy_str="";
										for(my $i=1;$i<= ($string_length-length($potential_previous_spacer));$i++)
											{
												$dummy_str=$dummy_str."-";
											}
										$potential_previous_spacer=$dummy_str.$potential_previous_spacer;
										
										#print "********** $potential_previous_spacer ********\n";
										#print "Skipping as: $cur_r_start-$string_length-1 <0 \n";
										#$skip_this_rec=1;last;
									}
								else{	
										if(($cur_r_start-$string_length-1)<length($species_seq))
											{						
												$potential_previous_spacer=substr($species_seq,($cur_r_start-$string_length-1),$string_length);		#	.$potential_previous_spacer;					
											}	
									}	
								
								#print "L: \$no_of_tries=$no_of_tries \$string_length=$string_length\n $potential_previous_spacer\n",length($potential_previous_spacer),"\n";
								#--------------------- now get the left flank ---------------------------------------------------
								#if(($cur_r_start-1-200)>0)
								#	{
								#		$left_flank=substr($species_seq,($cur_r_start-1-200),200);
								#	}
								#else{
								#		$left_flank=substr($species_seq,0,($cur_r_start-1));
								#	}
								#------------------------------------------------------------------------------------------------
								my($l_spacer,$l_spacer_pos);
								my $reference_repeat=$model_repeat;
								if($ea_dynamic_search==1 and $cur_r_seq!~/-/ and not $cur_no_of_insertions >0)
									{
										$reference_repeat=&change_dots_to_bases($cur_r_seq,$model_repeat);
									}
								($original_best_fitting_pattern,$best_fitting_pattern,$l_spacer,$l_spacer_pos,$insertion_base,$position)=&get_best_fitting_string_in_flanks($range,$accession,$allowed_percent_similarity,$potential_previous_spacer,$reference_repeat,"LEFT",1); # 1 will fill the gap with next bases
								
								
								#print "\nL. $original_best_fitting_pattern,$best_fitting_pattern,$l_spacer,$l_spacer_pos,$insertion_base,$position\nwith model_repeat:  $model_repeat\n";
								
								if(not defined $l_spacer)
									{
										##------ very important, don't block this or else the length of the string wont increase-----
										$best_fitting_pattern="-";
										next;
									}
								if($l_spacer eq "-"){$l_spacer="";}
								
								
								my $similarity_score;
								
								if(not $best_fitting_pattern)
									{
										$best_fitting_pattern=""; 										
										$similarity_score=0;
									}
								else{
										$similarity_score=&get_similarity_score($best_fitting_pattern,$model_repeat);
										#print "\$similarity_score=$similarity_score\n";
										#---- find out how many gaps are there ---
										my $no_gaps=$best_fitting_pattern=~s/-/-/g;
										if($no_gaps>1)
											{
												$similarity_score=$similarity_score -$no_gaps*2;
											}
										if($insertion_base)
											{
												my $no_of_ins= $insertion_base=~s/,/,/g;
												#print "\$insertion_base==$insertion_base\t $similarity_score=$similarity_score-$no_of_ins*2\t$similarity_score<(length($model_repeat)*$allowed_percent_similarity/100)\n";
												$similarity_score=$similarity_score-$no_of_ins*2;
											}
									}	
									
								#--------------------------------------------------------------------------------	
								#----- check best fitting pattern to have more than 60% similarity with $model repeat
								if($best_fitting_pattern eq "" or $similarity_score<(length($model_repeat)*$allowed_percent_similarity/100))
									{
										##------ very important, don't block this or else the length of the string wont increase-----
										$best_fitting_pattern="-";
										next;
									}
								if(length($l_spacer)>$max_gap_between_crisprs)#length($last_spacer_seq)*2.5)
									{		
										#print "Skipping...\n";
																		
										$skip_this_rec=1;#------ no need to push anything as the record is already there -----------------------
										last;
									}												
								
								
								
								#------------------------ now create the two lines which will be pushed
									
								#------------ increase the $extension_found -----------------------------------------
								$case_found=1;
								$no_of_tries=0;  #---reset the $no_of_tries
								$extension_found++;	
								$potential_previous_repeat=&change_bases_to_dots($best_fitting_pattern,$model_repeat);
								#------------------------------------------------------------------------------------
																
								
								#-------- now store the new rec ------------------------------------------------
								#my $gapless_p_r_seq=$potential_next_repeat; $gapless_p_r_seq=~s/-//g;
								my $gapless_p_r_seq=$potential_previous_repeat; $gapless_p_r_seq=~s/-//g;
								my $gapless_p_s_seq=$l_spacer; $gapless_p_s_seq=~s/-//g;						
								
								my $comma_removed_insertion_bases=$insertion_base;
								$comma_removed_insertion_bases=~s/,//g;
									
								my $previous_repeat_start=$cur_r_start-length($gapless_p_r_seq)-length($gapless_p_s_seq)-length($comma_removed_insertion_bases);	
									
								#my $potential_previous_spacer="";
								
								 
								if(not $l_spacer)
									{
										$l_spacer="";	
										$potential_previous_spacer="";									
										$pos=$r_start+$r_length;
										
										for(my $i=0;$i<length($last_spacer_seq);$i++)
											{
												$potential_previous_spacer=$potential_previous_spacer."-";
											}
									}
								else{
										$pos=$previous_repeat_start;		#-- deletion always happens just befor the current_start								
										$potential_previous_spacer=$l_spacer;
									}
									
								#$potential_previous_spacer=~s/-//g;						
								#--------------------------------------------------------------------------	
								
								
								
								#-------------------------------------------------------------------------  
								#-------- now get the actual base postions for the insertions -------------
								my @array_positions=split(',',$position);
								my $new_position="";
								foreach my $p(@array_positions)
									{
										my $p1=$previous_repeat_start+$p;
										$new_position=$new_position.",".$p1;
									}
								$new_position=~s/^,//;
								$position=$new_position;  
								my $comment2=""; 							
								if($insertion_base ne "")
									{
										$comment2="$insertion_base [$position]";	
									}
														
								my $new_rec_line1="$previous_repeat_start\t$potential_previous_repeat\t$potential_previous_spacer\t$comment2";
								#print "\$new_rec_line1=$new_rec_line1\n\n";
								
								unshift(@{$modified_array},$new_rec_line1);
								
								
								#---------------------------------- now shorten the flank----------------------------------------
								#my $gapless_best_fitting_pattern=$original_best_fitting_pattern;
								#   $gapless_best_fitting_pattern=~s/-//g;
								  
								  
								#my $chopping_pattern= $gapless_best_fitting_pattern.$l_spacer; chomp $chopping_pattern; $chopping_pattern=~s/\s+$//;
								
								#print "\$chopping_pattern=$chopping_pattern with \$potential_previous_spacer=$potential_previous_spacer\n";
								
								#$potential_previous_spacer=~s/$chopping_pattern$//;
								
								
								#print "\$potential_previous_spacer=$potential_previous_spacer\n\n";
								#$potential_previous_spacer=$left_flank;							
							} #--- end of while
						
						if($skip_this_rec==1){next;}
									
						
					}
				#------------ center contents ---------------	
				elsif($k1>4 and $k1<$#{$current_array}-1)
					{
						push(@{$modified_array},$current_line);	
						next;
					}
				#------------- right flank	
				else{
						push(@{$modified_array},$current_line);	
						
						#next;
						
						#print "Inside\n\n\n";
						#$s_seq=$right_flank;
						
						
						my @arr_previous_rec=split('\t',$$modified_array[$#{$modified_array}-1]);		#split('\t',$$current_array[$k1-1]);												
						
						#my $last_repeat_start=$arr_previous_rec[0];
						my $last_repeat_seq=$arr_previous_rec[1];
						#my $last_spacer_seq=$arr_previous_rec[2];  #--- remember, last_spacer_seq will always be 0

						my $last_spacer_seq;
						
						if(defined $arr_previous_rec[2])
							{
								$last_spacer_seq=$arr_previous_rec[2];  
							}
						else{
								$last_spacer_seq="";
							}
						
						#--- use while condition -----------------------------------------------------------
						my ($original_best_fitting_pattern,$best_fitting_pattern,$insertion_base,$position);
						
						$best_fitting_pattern="-";
						
						my $pos;
						#my $next_spacer_seq="";
						my $repeat_seq_before_flank=$r_seq;
						
						my $skip_this_rec=0;
						my $extension_found=0;
						my $potential_next_repeat="";
						my $potential_previous_spacer="";
						
						my $no_of_tries=0;
						while($best_fitting_pattern ne "")
							{
								
								$no_of_tries++;
								if($no_of_tries>5){last;}
								#-------- find the repeat_start position from the last record -----------
								my @arr_cur_rec=split('\t',$$modified_array[$#{$modified_array}]);	
								
								#print "\t@arr_cur_rec\n";
								my $cur_r_start=$arr_cur_rec[0];
								my $cur_r_seq=$arr_cur_rec[1];	my $gapless_cur_r_seq=$cur_r_seq; $gapless_cur_r_seq=~s/-//g;										
								my $cur_s_seq="";
								if($arr_cur_rec[2])
									{
										$cur_s_seq=$arr_cur_rec[2]; $cur_s_seq=~s/-//g;
									}
								
								#if($cur_s_seq !=0)
								#	{
								#		$skip_this_rec=1;
								#		last;
								#	}
								my $cur_comment=$arr_cur_rec[3];
								
								my $cur_no_of_insertions=0;
								my ($cur_insertion_bases,$cur_insertion_positions);
								if($arr_cur_rec[3])
									{
										#$comment=$tmp_array[3];
										my $c_comment=$cur_comment; $c_comment=~s/^\s+//;
										#if($c_comment=~/^Insertions/)
										#	{
										#		$c_comment=~ s/Insertions of //;
										#	}
										#els
										if($c_comment!~/^Del/)
											{										
												my @tmp_array=split(' ',$c_comment);
												
												#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
												$cur_insertion_bases=$tmp_array[0];
												$cur_insertion_positions=$tmp_array[1];
												#print "\$cur_insertion_bases=$cur_insertion_bases\n";
												
												$cur_insertion_bases=~s/,//g;
												$cur_no_of_insertions=length($cur_insertion_bases);
												#print "\$cur_no_of_insertions=$cur_no_of_insertions\n";
											}
									}							

								#--------------------------------------------------------------------------------------------------				
								my $string_length=0;
								
								if(length($last_spacer_seq)>15 and length($last_spacer_seq)<70)
									{
										$string_length=int((length($last_repeat_seq) + length($last_spacer_seq)*1.33*$no_of_tries));	#--- giving 33% extra
									}
								else{
										$string_length=int(length($last_repeat_seq)*2.5);	#--- giving 33% extra
									}	
									
									
								if($cur_r_start+length($cur_r_seq)+$cur_no_of_insertions-1+$string_length>length($species_seq))
									{
										
										$potential_previous_spacer="";
										
										#-------- make a dummy string using remaining bases and '-'s
										
										if(($cur_r_start+length($cur_r_seq)+$cur_no_of_insertions-1) < length($species_seq))
											{
												$potential_previous_spacer= substr($species_seq,($cur_r_start+length($gapless_cur_r_seq)+$cur_no_of_insertions-1),$string_length);
												#print "**** $r_start-1 : $previous_spacer_bases_to_chop\n\n";
											}	
										
										
										#-------- make a dummy string and take substr of it : currently just ending the process----------------- Example: NC_019776
										#print "\n\n############## here #################\n\n";
										my $dummy_str="";
										for(my $i=1;$i<= ($string_length-length($potential_previous_spacer));$i++)
											{
												$dummy_str=$dummy_str."-";
											}
										$potential_previous_spacer=$potential_previous_spacer.$dummy_str;
										
										
										#$skip_this_rec=1;last;
									}
								else{					
										$potential_previous_spacer=substr($species_seq,($cur_r_start+length($gapless_cur_r_seq)+$cur_no_of_insertions-1),$string_length);
									}
								#print "\t\$no_of_tries=$no_of_tries \t\$potential_previous_spacer=$potential_previous_spacer\n";
								
								
								#------- get the right_flank set as well -----------------------------------------------------------
								#if(($cur_r_start+length($cur_r_seq)-1+200)<length($species_seq))
								#	{
								#		$right_flank=substr($species_seq,($cur_r_start+length($cur_r_seq)+$cur_no_of_insertions-1),200);
								#	}
								#else{
								#		$right_flank=substr($species_seq,($cur_r_start+length($cur_r_seq)+$cur_no_of_insertions-1));
								#	}	
								#--------------------------------------------------------------------------------------------------
								
								my($l_spacer,$l_spacer_pos);
								
								my $reference_repeat=$model_repeat;
								if($ea_dynamic_search==1 and $cur_r_seq!~/-/ and not $cur_no_of_insertions >0)
									{
										$reference_repeat=&change_dots_to_bases($cur_r_seq,$model_repeat);
									}
									
								($original_best_fitting_pattern,$best_fitting_pattern,$l_spacer,$l_spacer_pos,$insertion_base,$position)=&get_best_fitting_string_in_flanks($range,$accession,$allowed_percent_similarity,$potential_previous_spacer,$reference_repeat,"RIGHT",1); # 1 will fill the gap with next bases
								#print "\nR: $no_of_tries : . Model_repeat:$reference_repeat\t $original_best_fitting_pattern, $best_fitting_pattern, $l_spacer, $l_spacer_pos, $insertion_base, $position\n\n";	
								
								if(not defined $l_spacer)
									{
										$best_fitting_pattern="-";
										next;
									}
								if($l_spacer eq "-"){$l_spacer="";}
								
								#print "OK1\n";
								my $similarity_score;
								
								if(not $best_fitting_pattern)
									{
										$best_fitting_pattern="";
										$similarity_score=0;
									}
								else{
										$similarity_score=&get_similarity_score($best_fitting_pattern,$model_repeat);
										
										#---- find out how many gaps are there ---
										my $no_gaps=$best_fitting_pattern=~s/-/-/g;
										if($no_gaps>1)
											{
												$similarity_score=$similarity_score -$no_gaps*2;
											}
										if($insertion_base)
											{
												my $no_of_ins= $insertion_base=~s/,/,/g;
												$similarity_score=$similarity_score-$no_of_ins*2;
											}
									}	
								#print "OK2\n";
								#--------------------------------------------------------------------------------	
								#----- check best fitting pattern to have more than 60% similarity with $model repeat
								if($best_fitting_pattern eq "" or $similarity_score<(length($model_repeat)*$allowed_percent_similarity/100))
									{
										##------ very important, don't block this -----
										$best_fitting_pattern="-";
										next;
									}
								if(length($l_spacer)>$max_gap_between_crisprs)#(length($last_spacer_seq)*1.5))
									{	
										#print "\$no_of_tries : Skipping...\n";									
										#print "\nR:$no_of_tries. Model_repeat:$reference_repeat\t $original_best_fitting_pattern, $best_fitting_pattern, $l_spacer, $l_spacer_pos, $insertion_base, $position\n\n";	
										
										$skip_this_rec=1;#------ no need to push anything as the record is already there -----------------------										
										#next;
										last;
									}												
								
								
								#print "OK4\n";
								#------------------------ now create the two lines which will be pushed
									
								#------------ increase the $extension_found -----------------------------------------
								$case_found=1;
								$no_of_tries=0;  #---reset the $no_of_tries
								$extension_found++;	
								$potential_next_repeat=&change_bases_to_dots($best_fitting_pattern,$model_repeat);
								#------------------------------------------------------------------------------------
								#my $pos;
								my $current_repeat=$cur_r_seq;		
								my $next_repeat=&change_bases_to_dots($best_fitting_pattern,$model_repeat);
								
								
								#-------- now store the new rec ------------------------------------------------
								#my $gapless_p_r_seq=$potential_next_repeat; $gapless_p_r_seq=~s/-//g;
								my $gapless_p_r_seq=$cur_r_seq; $gapless_p_r_seq=~s/-//g;
								my $gapless_p_s_seq=$l_spacer; $gapless_p_s_seq=~s/-//g;						
								
								
								my $comma_removed_insertion_bases=$insertion_base;
								$comma_removed_insertion_bases=~s/,//g;
									
								my $next_repeat_start=$cur_r_start+length($gapless_p_r_seq)+length($gapless_p_s_seq)+$cur_no_of_insertions;#+length($comma_removed_insertion_bases);	
								
								#print "$cur_r_start+length($gapless_p_r_seq)+length($gapless_p_s_seq)+$cur_no_of_insertions\n";
								#my $current_spacer=0;
								
								#print "OK\n";
								
								#--------------------------------------------------------------------------								
								my $comment1="";						
								if(length($l_spacer)<(length($last_spacer_seq)*80/100))
									{	
										if($cur_insertion_bases)
											{							
												$comment1=$cur_comment." Deletion [$next_repeat_start]";													
											}
										else{
												$comment1="Deletion [$next_repeat_start]";	
											}		
									}
								else{
										if($cur_insertion_bases)
											{							
												$comment1=$cur_comment;													
											}
										else{
												$comment1="";	
											}
									}	
								
								#------------ if l_spacer is blank, fill it with dashes -----------------
								if(not $l_spacer)
									{
										$l_spacer="";										
										$pos=$r_start+$r_length;
										
										#print "OK :",length($cur_s_seq),"\n";
										
										for(my $i=0;$i<length($last_spacer_seq);$i++)
											{											
												$l_spacer=$l_spacer."-";
											}
									}
								else{
										$pos=$next_repeat_start;		#-- deletion always happens just befor the next_start								
										$cur_s_seq=$l_spacer;
									}
									
								$cur_s_seq=~s/-//g;		
													
								#------- modify the previous record -------------------------------------
									
								my $modified_previous_rec="$cur_r_start\t$cur_r_seq\t$l_spacer\t$comment1";	
								$$modified_array[$#{$modified_array}]=$modified_previous_rec;
								
								
								
								#-------------------------------------------------------------------------  
								#-------- now get the actual base postions for the insertions -------------
								my @array_positions=split(',',$position);
								my $new_position="";
								foreach my $p(@array_positions)
									{
										my $p1=$next_repeat_start+$p;
										$new_position=$new_position.",".$p1;
									}
								$new_position=~s/^,//;
								$position=$new_position;  
								my $comment2=""; 							
								if($insertion_base ne "")
									{
										$comment2="$insertion_base [$position]";	
									}
														
								my $new_rec_line1="$next_repeat_start\t$potential_next_repeat\t0\t$comment2";
								#print "\$new_rec_line1=$new_rec_line1\n\n";
								
								push(@{$modified_array},$new_rec_line1);
								
								#print "OK5\n";
								#---------------------------------- now shorten the flank----------------------------------------
								#my $gapless_best_fitting_pattern=$original_best_fitting_pattern;
								#   $gapless_best_fitting_pattern=~s/-//g;
								   
								#my $chopping_pattern= $l_spacer.$gapless_best_fitting_pattern; 
								
								#$right_flank=~s/^$chopping_pattern//;
								#$s_seq=$right_flank;
								
								#exit;							
							} #--- end of while
						
						if($skip_this_rec==1){next;}
									
					
					}
													
			}
		return($case_found);		
	}



sub get_best_fitting_string_in_flanks()
	{
		my($range,$accession,$allowed_percent_similarity,$bases,$model_repeat_bases,$side,$fill_trailing_gaps)=@_;
		
		if(not defined $bases or $bases eq ""){return("","","","","","");}
		
		$bases=~s/^-+//;
		$bases=~s/-+$//;
		
		#print "\nInside sub: $accession,$allowed_percent_similarity,$bases,$model_repeat_bases,$side\n"; #return;
		
		my $original_best_fitting_pattern="";
		my $best_fitting_pattern="";
		
		#my $model_repeat_bases="";
		my $insertion_bases="";
		my $insertion_bases_positions="";
		my $last_spacer="";
		my $last_spacer_pos="";
					
		my @array_bases;
		my @array_model_repeat_bases;
		
		
		
		@array_bases=split('',$bases);	
		@array_model_repeat_bases=split('',$model_repeat_bases);
		
		
		#my $letters="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
		#my @arr_letters=split('',$letters);
		
		
		
		
		#---- now get pairwise alignments using water
		#if($no_of_gaps>1)
		#	{

						
				#for(my $i=$no_of_gaps;$i>0;$i--)
				#	{
					my $top_line="";
					my $bottom_line="--";
					
					my $n_o_g=0;
					my $n_o_insertion=0;
					my $difference_in_start_position=0;
					my $gap_increase=0;
					
					my $best_top_line="";
					my $best_bottom_line="";
					
					my $no_of_tries=0;
					my %hash_of_alignments;
					#system("echo 'Going into whhile1 with $bases,$model_repeat_bases' >>log1.txt");
					while($no_of_tries<5)
					{
						$no_of_tries++;
						
						my $gapopen=4.5+$gap_increase;
						my $gapextend=$gapopen-0.05; 						
						if($gapextend>10){$gapextend=9.5;}
						
						#system("echo 'Inside whhile1 with $bases,$model_repeat_bases' >>log1.txt");
						
						my $outfile=&run_water($range,$accession,$bases,$model_repeat_bases,$gapopen,$gapextend);
						#---- now open the output file and get the alignment, store the alignment in a compiled file
						
						
						if(-e "$tmp_dir\/$outfile")
							{
								open(RD1,"$tmp_dir\/$outfile") or print "";#ERROR in __LINE__ $! - can't find $tmp_dir\/$outfile\n<br>";
								flock(RD1,2);
								my @arr_rd=<RD1>;
								close(RD1);
								unlink("$tmp_dir\/$outfile");
								
								
								my $tl_count=0;
								my $bl_count=0;
								foreach my $line(@arr_rd)
									{
										#if($line=~/TACGGTACATTAGGAAC/)
										#{
										#print $line,"\n";
										#}
										chomp $line;$line=~s/\r//g;
										if(not $line or $line=~/#/){next;}
										elsif($line=~/\d+/)
											{
												#print $line,"\n";
												$line=~s/^\s+//;$line=~s/\s+/\t/g;
												my($start,$seq,$stop)=split('\t',$line);
												if($tl_count==0){$top_line=$line;$tl_count++;}
												elsif($bl_count==0){$bottom_line=$line;$bl_count++;}
											}	
										
										#if($line=~/score/i)
										#	{
										#		print "$line\n";
										#	}
									}
								
								#print "\n";	

								
						
							}
						#else{next;}
						
						#print "$top_line\n$bottom_line\n\n<br>";
						
						if($top_line eq "" or $bottom_line eq ""){next;}
						
						my($top_start,$top_seq,$top_stop)=split('\t',$top_line);
						my($bottom_start,$bottom_seq,$bottom_stop)=split('\t',$bottom_line);
						
						if(not defined $top_seq or not defined $bottom_seq){next;}
						
						#--------------------------------------------------------------------
						
						$n_o_g=$bottom_seq=~ s/-/-/g;
						if(not $n_o_g){$n_o_g=0;}
						
						$n_o_insertion=$top_seq=~ s/-/-/g;
						if(not $n_o_insertion){$n_o_insertion=0;}
						
						
						$difference_in_start_position=$top_start-$bottom_start;
						
						
						#------ now calculate the score and store in an array
						my $similarity_score=&get_similarity_score($top_seq,$bottom_seq);
						$similarity_score=$similarity_score;#-$n_o_g-$n_o_insertion**2;
						
						$hash_of_alignments{$similarity_score}{$n_o_g}=$top_line.",".$bottom_line;
						
						#print "\$similarity_score=$similarity_score\t\$n_o_g=$n_o_g\n\n";
						
						
						
						#if($n_o_g<=0 )      # max one insertion is allowed, ($top_start-$bottom_start)<=1 will achieve that
						#	{
						#		$best_top_line=$top_line;
						#		$best_bottom_line=$bottom_line;
						#		last;
						#	}
						#els
						if($no_of_tries>=5)
							{
								
								#---- now get the best scorring string---------
								my $best_rec_index=0;
								my $best_score=0;
								my $best_rec_line="";
								my $best_rec_nog=0;
								foreach my $s_score(sort{$b<=>$a}keys %hash_of_alignments)
									{
										foreach my $gaps(sort{$a<=>$b}keys %{$hash_of_alignments{$s_score}})
											{
												#print "$best_rec : $s_score\t$gaps\t$hash_of_alignments{$s_score}{$gaps}\n";
												
												
												$best_score=$s_score;
												$best_rec_nog=$gaps;
												$best_rec_line=$hash_of_alignments{$s_score}{$gaps};
												
												$best_rec_index++;
												last;
											}
										if($best_rec_index>0){last;}	
									}
								#----------------------------------------------
								
								if($best_rec_index>0 and $best_score>(length($model_repeat_bases)*($allowed_percent_similarity)/100)) 
									{
										($best_top_line,$best_bottom_line)=split(',',$best_rec_line);
										
									}
								else{	
										#print "skipping from sub with: $best_rec_line\n";
										return("","","","","","");
										#-------------------------------------------------------
										#$best_top_line="0\t$bases\t$no_of_gaps";
										#$best_bottom_line="0\t$model_repeat_bases\t$no_of_gaps";
										$best_bottom_line="-";
									}
									
								#---------------------------------------------------------------	
								$n_o_g=0;
								$difference_in_start_position=0;
								last;
								
							}	
						#$no_of_tries++;
						$gap_increase++;
					}	
						
					#system("echo 'Outside whhile1 with $bases,$model_repeat_bases' >>log1.txt");	
					#print "\n\nBEST TOP LINE: $best_top_line\nBest Bottom line: $best_bottom_line\n";	
					
					

					
					#---- read this first
					#-- not working with ($top_start-$bottom_start)<=1 ; may need to individually delete each base at a time, and check the alignment
					#---- read above line
					#--------------------------------------------------
					
					
					#----- now complete check the alignment and add bases which is out side the alignment
					
					# Case 1: no gaps in the bottom line
					# case 2: one gap in the bottom line	
					# case 3: no alignment found	
					
					if(not defined $best_top_line or not defined $best_bottom_line){return("","","","","","");}
					
					my($top_start,$top_seq,$top_stop)=split('\t',$best_top_line);
					if(not defined $top_start or not defined $top_stop){return("","","","","","");}
					
					$top_start=$top_start-1;
					$top_stop=$top_stop-1;
					
					my($bottom_start,$bottom_seq,$bottom_stop)=split('\t',$best_bottom_line);
					if(not defined $bottom_start or not defined $bottom_stop){return("","","","","","");}
					$bottom_start=$bottom_start-1;
					$bottom_stop=$bottom_stop-1;
					
					if(not defined $top_seq or not defined $bottom_seq){return("","","","","","");}
					
					#print "\$difference_in_start_position=$difference_in_start_position\n";
					
					
					#--- read this first:
					#--- make the string in such way that there will be Gaps either end and finally just use the model string's length to use substr
					#-- for LEFT and RIGHT the conditions will vary
					# case 1: the matching bases can be exactly in the middle and there can be unmatched bases in both side
					# case 2: matching bases are on the left
					# case 3: matching bases are on the right
					# case 4: Insertion base(s) are the bases which have  higher position index than the length of the model repeat string after gap insertion, and if a natural insertion would have occured in water alignment
					
					
					
					#------------ now handle the inserted bases --------------------------------
					my @arr_top_seq=split('',$top_seq);
					my @arr_bottom_seq=split('',$bottom_seq);
					
					
					for(my $i=0;$i<=$#arr_bottom_seq;$i++)
						{
							if($arr_bottom_seq[$i] eq "-")
								{
									if(defined $arr_bottom_seq[$i-1] and $arr_bottom_seq[$i-1] eq "-")
										{
											$insertion_bases=$insertion_bases.$arr_top_seq[$i];
											#$insertion_bases_positions=$insertion_bases_positions.$i;
										}
									else{	
											$insertion_bases=$insertion_bases.",".$arr_top_seq[$i];
											$insertion_bases_positions=$insertion_bases_positions.",".$i;
										}
									$arr_top_seq[$i]="";
								}
						}
					$insertion_bases=~s/^,//;	
					$insertion_bases_positions=~s/^,//;
					my $insertion_removed_bfp=join("",@arr_top_seq);
					
					#print "$insertion_removed_bfp\t$insertion_bases\t$insertion_bases_positions\n\n";
					#---------------------------------------------------------------------------
					
					#my @arr_bfp=split("",$top_seq);
					my @arr_bfp=split("",$insertion_removed_bfp);
					my @arr_bfp_original=split("",$top_seq);
					
					my $index_1=1;
					my @tmp_arr1;
					for(my $k1=$bottom_start-1;$k1>=0;$k1--)
						{
							if(($top_start-$index_1)>=0)
								{
									unshift(@tmp_arr1,$array_bases[$top_start-$index_1]);
								}
							else{
									unshift(@tmp_arr1,"-");
								}	
							
							$index_1++;
						}
									
									
					my $index_2=1;
					my @tmp_arr2;
									
					#print "\$bottom_stop=$bottom_stop \$#array_model_repeat_bases=$#array_model_repeat_bases\n";								
					
					
					
					#---- only fill remaining gaps when extending in the flanks 
					if($fill_trailing_gaps==1)
						{			
							for(my $k2=$bottom_stop+1;$k2<=$#array_model_repeat_bases;$k2++)
								{
															#print "Inside \$array_bases[$top_stop+$index_2]=$array_bases[$top_stop+$index_2]\n";
									if(($top_stop+$index_2)<=$#array_bases and defined $array_bases[$top_stop+$index_2])
										{
											push(@tmp_arr2,$array_bases[$top_stop+$index_2]);
										}			
									else{
											push(@tmp_arr2,"-");
										}
									$index_2++;
								}
						
						}
					#------- get best_fitting_pattern without the insertions	
					#if(not defined $tmp_arr1[0]){return("","","","","","");} #---- never enable this, as then the arrays will not be extended 	in the flanks.
					if(defined $tmp_arr1[0])
						{
							unshift(@arr_bfp,join("",@tmp_arr1));
						}				
					if(defined $tmp_arr2[0])
						{
							push(@arr_bfp,join("",@tmp_arr2));									
						}	
					$best_fitting_pattern=join("",@arr_bfp);	
					#print "\$best_fitting_pattern=$best_fitting_pattern\n";	
					
					#--- now get another BFP with the insertion, this requires to chop the beginning or ends of the flanks
					if(defined $tmp_arr1[0])
						{
							unshift(@arr_bfp_original,join("",@tmp_arr1));
						}
					if(defined $tmp_arr2[0])
						{		
							push(@arr_bfp_original,join("",@tmp_arr2));	
						}									
					$original_best_fitting_pattern=join("",@arr_bfp_original);	
					
					#print "\$original_best_fitting_pattern=$original_best_fitting_pattern\n";	
					
					
					if($side eq "LEFT") #  all operations including base deletion should start on the left side
						{	
									#print join("",@array_bases),"\n";
									#---------- now get the insertion bases and position
									my $gap_less_bfp=$original_best_fitting_pattern;
									$gap_less_bfp=~s/-//g;
									
									
									if($bases=~s/$gap_less_bfp$//)
										{
											#$insertion_base="";
											#$insertion_base_position=0;
											
											$last_spacer="";
											$last_spacer_pos=0;
											
											
										}
									else{
											#print "inside\n";
											my @tmp_arr3=split($gap_less_bfp,$bases);
											
											#$insertion_base=$tmp_arr3[1];
											#$insertion_base_position=length($gap_less_bfp);
											
											$last_spacer=$tmp_arr3[$#tmp_arr3]; 			#--- the last record holds the last spacer
											$last_spacer_pos=length($gap_less_bfp);
										}
									#$insertion_base=substr($bases,0,abs($#array_bases-$#array_model_repeat_bases));
									#$insertion_base_position=0;
									
									
									
									#print "\nA:	\$best_fitting_pattern=$best_fitting_pattern \t\$last_spacer=$last_spacer\n";	
								
						
								
						}
					
					elsif($side eq "RIGHT") #  all operations including base deletion should start on the right side
						{				
									
									#print join("",@array_bases),"\n";
									#---------- now get the insertion bases and position
									my $gap_less_bfp=$original_best_fitting_pattern;
									$gap_less_bfp=~s/-//g;
									
									
									if($bases=~s/^$gap_less_bfp//) #--- for the cases where the whole spacer got deleted
										{
											$last_spacer="";
											$last_spacer_pos=0;
											
											#$insertion_base="";
											#$insertion_base_position=0;
										}
									else{						#--- for the cases where the part of spacer got deleted
											
											my @tmp_arr3=split($gap_less_bfp,$bases);
											
											#$insertion_base=$tmp_arr3[0]; 
											#$insertion_base_position=0;
											
											$last_spacer=$tmp_arr3[0]; 			#--- the first record holds the last spacer
											$last_spacer_pos=0;
										}
												
									
									
									#print "\nB:	\$best_fitting_pattern=$best_fitting_pattern with $insertion_base\t found at $insertion_base_position\n";	
								
						}
					
	#		}
	#	else{
	#			return("","","","","","");
				#if($side eq "LEFT")
				#	{	
				#		$best_fitting_pattern=$array_bases[$#array_bases];							
				#	}
				#elsif($side eq "RIGHT")
				#	{
				#		$best_fitting_pattern=$array_bases[0];
				#	}
				
	#		}	
		if($last_spacer eq ""){$last_spacer="-";}
		#print "$original_best_fitting_pattern,$best_fitting_pattern,$last_spacer,$last_spacer_pos,$insertion_bases,$insertion_bases_positions\n\n";
		return ($original_best_fitting_pattern,$best_fitting_pattern,$last_spacer,$last_spacer_pos,$insertion_bases,$insertion_bases_positions);
	}







sub fix_arrays_with_gaps()
	{
		
		my($range,$accession,$model_repeat,$avg_spacer_length,$current_array,$modified_array)=@_;
		
		#print "Going to fix repeats with gap(s)  for $accession :\n";
		#system("echo 'Going to fix repeats with gap(s)  for $accession :' >log.log");
		
		my $case_found=0;
		my $left_flank="";
		my $right_flank="";
		
		#print "Original \$left_flank=$left_flank\n";
		#print "Original \$right_flank=$right_flank\n";
		
		
		
		#---------- now open the sequence file and get the sequence string -------------
		open(SEQ,"$tmp_dir\/$accession\.fna") or print "$!";
		my @arr_seq=<SEQ>;
		close(SEQ);
		my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;
		
		
		
		my $coord_diff=0;
		
		
		
		my $first_occurrence_of_gap=0;
		
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
					{
						#print "@{$current_array-[0]}\n";
						my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
							
						#print "$current_line\n<br>";
						#system("echo 'Line:$k1-$current_line' >>log.log");
							
						my @tmp_array=split('\t',$current_line);
						my($r_start,$r_length,$s_length,$r_seq,$s_seq);
						my $comment="";	
						my $insertion_base_and_position="";
						
						
						$r_start=$tmp_array[0];
						

						
						$r_seq=$tmp_array[1]; $r_seq=~s/\s+//g; my $r_seq1=$r_seq; $r_seq1=~s/-//g;
						if(not defined $tmp_array[2])
							{
								$s_seq=""; 
							}
						else{	
								$s_seq=$tmp_array[2]; 
							}
						$s_seq=~s/\s+//g; 
						
						my $s_seq1=$s_seq; $s_seq1=~s/-//g;
								
						my $existing_insertion_bases="";
						my $existing_insertion_positions="";	
						if($tmp_array[3])
							{
								$comment=$tmp_array[3];
									
								#print "\n\n\$comment=$comment\n";
								if($tmp_array[3]!~/^\s{0,1}Del/)
									{		
										($existing_insertion_bases,$existing_insertion_positions)=split(' ',$comment);
										if($existing_insertion_positions)
											{
												$existing_insertion_positions=~s/\[//g;
												$existing_insertion_positions=~s/\]//g;
											}
									}		
									
								#print "\n\n\$comment=$comment\t\$existing_insertion_positions=$existing_insertion_positions\t\$existing_insertion_positions=$existing_insertion_positions\n";	
							}
								
						$r_length=length($r_seq1);
						$s_length=length($s_seq1);
						
						
						#----------------------- pre record ------------------------------------		
						my $pre_comment="";
						my $pre_spacer="";
						if($k1>4)
							{
								# ----get the previous rec_line ----------------------------------
								my @arr_pre_rec=split('\t',$$modified_array[$#{$modified_array}]);
								if(defined $arr_pre_rec[2])
									{
										$pre_spacer=$arr_pre_rec[2]; $pre_spacer=~s/-//g;
									}
								else{
										$pre_spacer="";
									}		
								
								if(defined $arr_pre_rec[3])
									{
										$pre_comment=$arr_pre_rec[3];
									}
							}
						 my $pre_spacer1=$pre_spacer; $pre_spacer1=~s/-//g;	
						#--------- determine from where to take the bases to fill the gaps---------
						
						my $take_seq_from="";
						my $gap_string="";
						if($r_seq=~/-/)
							{								
								$case_found=1;
								my ($best_fitting_pattern,$insertion_base,$position);	
										
								
								# case 1: when there is gaps in the very first repeat (on the left side), take bases from the left flank, if gaps on the right side, then take from spacer seq
								#-- no of bases to take is equal or less than the (no_of_gaps-1) present in the current repeat
								#--- so find the gaps position
								
								if($r_seq=~ /^-/)
									{
										$take_seq_from="LEFT";
										
										
										my $no_of_leading_gaps;
										if($r_seq=~ /^(\-+)/)
											{
												$gap_string=$1;
												
												#print "\$1=$1\n";
												$no_of_leading_gaps =length($1);
											}
											
											
										$best_fitting_pattern="-";
										
										#------------------------------------------------------------------------------
										
										#---- now take bases from the right side end of the flank/previous spacer
										#if($#{$modified_array}<0)
										if($k1==4)
											{
												
												
												
												my $factor_1=$no_of_leading_gaps+1;
												my $model_repeat_bases=substr($model_repeat,0,$no_of_leading_gaps); 
												
												my $no_of_tries=0;
												while($best_fitting_pattern=~/^-/)
													{
														
														$no_of_tries++;
														if($no_of_tries>10){last;}
														#---- left flanks are normally correct, the right flank is messed up when there are gaps in the CRISPR array
														
														#$left_flank=substr($species_seq,($cur_r_start-1-200),200);
														#$right_flank=substr($species_seq,($cur_r_start+length($cur_r_seq)+$cur_no_of_insertions-1),200);
														
														#my $bases=substr($left_flank,-($no_of_leading_gaps+$factor_1));		
														
														my $bases="";
														
														if(($r_start-1-$factor_1)>0)
															{
																$bases=substr($species_seq,($r_start-1-$factor_1),$factor_1);											
															}											
														
														#print "taken $no_of_leading_gaps ($bases) from the left flank : \n";										
															
														
														($best_fitting_pattern,$insertion_base,$position)=&get_best_fitting_string($range,$accession,$no_of_leading_gaps,$bases,$model_repeat_bases,"LEFT");
														
														#print "$best_fitting_pattern,$insertion_base,$position\n\n";							
												
														if($best_fitting_pattern eq ""){last;}
												
														
														my $no_of_l_gaps=0;
														if($best_fitting_pattern=~ /^(\-+)/)
															{
																#my $no_of_l_gaps=$1;
																
																#print "\$1=$1\n";
																$no_of_l_gaps =length($1);
															}
														#if($no_of_l_gaps>0 and length($insertion_base)>$no_of_l_gaps)
														#	{
														#		$best_fitting_pattern=substr($species_seq,$r_start-1-$no_of_leading_gaps,$no_of_leading_gaps);
														#		$insertion_base="";
														#		$position="";
														#		#print "$best_fitting_pattern,$insertion_base,$position\n\n";
														#		
														#		last;
														#	}
														
														#if(length($s_seq)>$avg_spacer_length)
														#	{
														#		$factor_1=$factor_1+(length($s_seq)-$avg_spacer_length);
														#	}
														#else{
															if($no_of_l_gaps)
																{
																	$factor_1=$factor_1+$no_of_l_gaps;
																}
															else{
																	$factor_1=$factor_1+1;
																}		
														#	}	
														if(($factor_1>(length($gap_string)*4)) or $factor_1>$s_length)
															{
																$best_fitting_pattern=substr($species_seq,$r_start-1-$no_of_leading_gaps,$no_of_leading_gaps);
																$insertion_base="";
																$position="";
																last;
															}
														
													}
												
												if($best_fitting_pattern eq ""){$best_fitting_pattern=$gap_string;}
												
												#-------------- now change identical bases with dots in $best_fitting_pattern
												$best_fitting_pattern=&change_bases_to_dots($best_fitting_pattern,$model_repeat_bases);
												
												my $gapless_best_fitting_pattern=$best_fitting_pattern;$gapless_best_fitting_pattern=~ s/-//g;
												my $no_of_bases_to_chop=length($gapless_best_fitting_pattern)+length($insertion_base);
												
												
												#my $chopping_pattern=substr($left_flank,-$no_of_bases_to_chop);
													
													
													
													#print "\$best_fitting_pattern=$best_fitting_pattern\n";
													#my $gapless_best_fitting_pattern=$best_fitting_pattern;$gapless_best_fitting_pattern=~ s/-//g;
													
													#--- now modify the previous repeat
													
													#print "\nOriginal \$left_flank : $left_flank\n\n";
												
												#$left_flank=~ s/$chopping_pattern$//;
												
												
												#---- modify the position  $r_start
												
												$r_start=$r_start-$no_of_bases_to_chop;		#-length($chopping_pattern);
												
												if($insertion_base ne "")
													{
														my $pos=$position;
														$position=$r_start+$position;
														
														if($existing_insertion_bases)
															{
																$insertion_base=$insertion_base.",".$existing_insertion_bases;
																$position=$position.",".$existing_insertion_positions;
															}
														$comment="$insertion_base [$position]";
													}
												elsif($existing_insertion_bases)
													{
														
														$comment="$existing_insertion_bases [$existing_insertion_positions]";
													}		
												
												$r_seq=~ s/^-+/$best_fitting_pattern/;
												
												#print "\nModified \$left_flank : $left_flank\n\n";
												
											}
										#elsif(length($pre_spacer1)>$avg_spacer_length*0.90)
										else{
												
												#system("echo 'Line:$k1' >log.log");
											
												#-------- check for base(s) deletion and skip if deletion is true -------------										
												
												
												my @arr_previous_rec=split('\t',$$current_array[$k1-1]);												
												my $last_spacer_seq=$arr_previous_rec[2];
												if(not defined $last_spacer_seq){$last_spacer_seq="";}
												#----- instead of using whole $$modified_array[$#{$modified_array}] line as a string, split and get the previous spacer
																								
												my $factor_1=$no_of_leading_gaps+1;
												my $model_repeat_bases=substr($model_repeat,0,$no_of_leading_gaps);												
												
												
												#if( (($r_length-$no_of_leading_gaps)<($r_length*3/4)) and (length($last_spacer_seq)<$avg_spacer_length))
												
												if($last_spacer_seq=~/$model_repeat_bases$/)
													{
														#print "$model_repeat_bases found in the previous spacer...\n";
													}
												elsif( ($last_spacer_seq!~/[ATGC]/) or ((($no_of_leading_gaps>=$r_length*10/100)) and (length($last_spacer_seq)<($avg_spacer_length-$avg_spacer_length*10/100))))
													{									
														#print "Skipping.. $current_line\n";
														my $new_rec_line="$r_start\t$r_seq\t$s_seq\t$comment";
														push(@{$modified_array},$new_rec_line);
														
														my $tmp_pos=$r_start;
														$comment="Deletion [$tmp_pos]";
														if(defined $tmp_arr_previous_rec1[0] and $tmp_arr_previous_rec1[0]=~/\S+/)
														{
														my @tmp_arr_previous_rec1=split('\t',$$modified_array[$#{$modified_array}-1]);
														$tmp_arr_previous_rec1[3]=$comment;
														
														my $old_rec_line=join(';',@tmp_arr_previous_rec1);
														$old_rec_line=~s/;/\t/g;
														$$modified_array[$#{$modified_array}-1]=$old_rec_line;
														}
														
														next;
													}
											

												
												#print "\n\n\n$model_repeat_bases\t$model_repeat\n\n\n"; 
												#print "FACTOR: $factor_1<br>";
												
												my $no_of_tries=0;
												while($best_fitting_pattern=~/^-/)
													{														
														
														$no_of_tries++;
														if($no_of_tries>10){last;}
														#print "\$last_spacer_seq=$last_spacer_seq\n\$no_of_leading_gaps=$no_of_leading_gaps\n";
														
														my $bases=substr($last_spacer_seq,-$factor_1); 
														
														
														
														
														#print "\ntaken $no_of_leading_gaps ( $bases ) from the previous spacer\n";
														
														($best_fitting_pattern,$insertion_base,$position)=&get_best_fitting_string($range,$accession,$no_of_leading_gaps,$bases,$model_repeat_bases,"LEFT");
														
														#print "$best_fitting_pattern,$insertion_base,$position\n\n";
														
														if($best_fitting_pattern eq ""){last;}
														
														#system("echo 'taken $no_of_leading_gaps ( $bases ) from the previous spacer $best_fitting_pattern,$insertion_base,$position' >>log.log");
														my $no_of_l_gaps=0;
														if($best_fitting_pattern=~ /^(\-+)/)
															{
																#my $no_of_l_gaps=$1;
																
																#print "\$1=$1\n";
																$no_of_l_gaps =length($1);
															}
														#if($no_of_l_gaps>0 and length($insertion_base)>$no_of_l_gaps)
														#	{
														#		$best_fitting_pattern=substr($last_spacer_seq,-$no_of_leading_gaps);
														#		$insertion_base="";
														#		undef $position;
														#		#print "$best_fitting_pattern,$insertion_base,$position\n\n";
														#		
														#		last;
														#	}
														
														
														#$factor_1=$factor_1+$no_of_l_gaps;
														
														if(length($s_seq)>$avg_spacer_length and ($#{$current_array}-4)>3) #------ just to avoid the smaller (2 spacer arrays)
															{
																$factor_1=$factor_1+(length($s_seq)-$avg_spacer_length);
															}
														else{
																if($no_of_l_gaps)
																	{
																		$factor_1=$factor_1+$no_of_l_gaps;
																	}
																else{
																		$factor_1=$factor_1+1;
																	}		
															}
														my $len=length($gap_string);	
														if(($factor_1>(length($gap_string)*4)) or $factor_1>$s_length )
															{
																#system("echo 'Final $factor_1 $best_fitting_pattern,$insertion_base,$position' >>log.log");
																$best_fitting_pattern=substr($last_spacer_seq,-$no_of_leading_gaps);
																$insertion_base="";
																undef $position;
																last;
															}
														
													}
												#print "$best_fitting_pattern,$insertion_base,$position\n\n";	
												if($best_fitting_pattern eq ""){$best_fitting_pattern=$gap_string;}
												
												#-------------- now change identical bases with dots in $best_fitting_pattern
												$best_fitting_pattern=&change_bases_to_dots($best_fitting_pattern,$model_repeat_bases);
												
												
																					
												my $gapless_best_fitting_pattern=$best_fitting_pattern;$gapless_best_fitting_pattern=~ s/-//g;
												my $no_of_bases_to_chop=length($gapless_best_fitting_pattern)+length($insertion_base);
												
												
												if($insertion_base ne "")
													{
														$position=$r_start-$no_of_bases_to_chop+$position;
													}
												#chomp $$modified_array[$#{$modified_array}]; $$modified_array[$#{$modified_array}]=~s/[\r|\n]$//;
												
												my $chopping_pattern=substr($last_spacer_seq,-$no_of_bases_to_chop);
												#print "\$no_of_bases_to_chop=$no_of_bases_to_chop\t\$chopping_pattern=$chopping_pattern from $last_spacer_seq\n";
											
												
												#--- now modify the previous spacer and update the record
												$last_spacer_seq=~ s/$chopping_pattern$//;
												
												$$modified_array[$#{$modified_array}]="$arr_previous_rec[0]\t$arr_previous_rec[1]\t$last_spacer_seq\t$insertion_base_and_position";
												$$current_array[$k1-1]="$arr_previous_rec[0]\t$arr_previous_rec[1]\t$last_spacer_seq\t$insertion_base_and_position";												
												

												
												
												
												#----now rectify the current record
												$r_start=$r_start-length($chopping_pattern);
												
												if($insertion_base)
													{
														#$position=$r_start+$position;
														if($existing_insertion_bases)
															{
																$insertion_base=$insertion_base.",".$existing_insertion_bases;
																$position=$position.",".$existing_insertion_positions;
															}
															
														$comment="$insertion_base [$position]";
														#print "<br>$comment<br>";
														
													}
												elsif($existing_insertion_bases and $existing_insertion_positions)
													{
														
														$comment="$existing_insertion_bases [$existing_insertion_positions]";
														
													}		
												
												
												
												$r_seq=~ s/^-+/$best_fitting_pattern/;
												
												#print "modified previous spacer : $$modified_array[$#{$modified_array}]\n";
												
												}
												
											
										
										#---- now replace the leading gaps with bases from the flank/previous spacer											
									}

								if($r_seq=~ /-$/)
									{
										$take_seq_from="RIGHT";
										
										my $no_of_trailing_gaps;
										if($r_seq=~ /(\-+)$/)
											{
												$gap_string=$1;
												#print "\$1=$1\n";
												$no_of_trailing_gaps =length($1);
											}
										
										
										$best_fitting_pattern="-";
											
										#my ($best_fitting_pattern,$insertion_base,$position);	
										#	$best_fitting_pattern="";	
											
										#---------------------------
										#---- now take bases from the right side end of the flank/previous spacer
										if($k1==($#{$current_array}-1))
											{
																								
												#--- now get the right flank after checking the corrds
												
												my $tmp_r_seq=$r_seq;
												$tmp_r_seq=~s/-//g;
												
												#my $last_record=$$modified_array[$#{$modified_array}];
												my $last_record=$$current_array[$k1-1];		
												#print "\$last_record=$last_record\n";	
																	
												my($last_rec_start,$last_repeat,$last_spacer,$memo_field)=split('\t',$last_record);
												
												$last_repeat=~s/-//g;
												$last_spacer=~s/-//g;
												
												$r_start=$last_rec_start+length($last_repeat)+length($last_spacer);
												
												my $r_flank_start=$r_start+length($tmp_r_seq);
												
												#print "\nOriginal Right Flank: $right_flank\n";
												
												if(($r_flank_start-1)<length($species_seq))
												{
												if(length($species_seq)> ($r_flank_start-1+100))
													{
														$right_flank=substr($species_seq,$r_flank_start-1,100);
													}
												else{
														$right_flank=substr($species_seq,$r_flank_start-1,(length($species_seq)-$r_flank_start));
													}
												}	
												#print "\nCorrected Right Flank: $right_flank\n<br>";
												#------------------------------------------------------------------------------------------------------
												
												
												
												my $factor_1=$no_of_trailing_gaps+1;
												my $model_repeat_bases=substr($model_repeat,-$no_of_trailing_gaps); # this one will be one base short from the $bases to allow 1 base insertion
												
												my $no_of_tries=0;		
												while($best_fitting_pattern=~/-$/)
													{														
														
														$no_of_tries++;
														if($no_of_tries>10){last;}
								
														my $bases=substr($right_flank,0,$factor_1);# this one will be one base longer that the $no_of_trailing_gaps 
														
														
														#print "taken $no_of_gaps ($bases ) from the left flank : \n";
														
														($best_fitting_pattern,$insertion_base,$position)=&get_best_fitting_string($range,$accession,$no_of_trailing_gaps,$bases,$model_repeat_bases,"RIGHT");
														
														#print "$best_fitting_pattern,$insertion_base,$position\n\n";
														if($best_fitting_pattern eq ""){last;}
												
														my $no_of_t_gaps=0;
														if($best_fitting_pattern=~ /(\-+)$/)
															{
																#my $no_of_l_gaps=$1;
																
																#print "\$1=$1\n";
																$no_of_t_gaps =length($1);
															}
															
														#if($no_of_t_gaps>0 and length($insertion_base)>$no_of_t_gaps)
														#	{
														#		$best_fitting_pattern=substr($right_flank,0,$no_of_trailing_gaps);
														#		$insertion_base="";
														#		undef $position;
														#		#print "$best_fitting_pattern,$insertion_base,$position\n\n";
														#		
														#		last;
														#	}	
														#$factor_1=$factor_1+$no_of_t_gaps;
														#if(length($s_seq)>$avg_spacer_length)
														#	{
														#		$factor_1=$factor_1+(length($s_seq)-$avg_spacer_length);
														#	}
														#else{
																if($no_of_t_gaps)
																	{
																		$factor_1=$factor_1+$no_of_t_gaps;
																	}
																else{
																		$factor_1=$factor_1+1;
																	}		
														#	}
														if(($factor_1>(length($gap_string)*4)) or $factor_1>$s_length)
															{
																$best_fitting_pattern=substr($right_flank,0,$no_of_trailing_gaps);
																$insertion_base="";
																undef $position;
																last;
															}
														
													}
												
												
												if($best_fitting_pattern eq ""){$best_fitting_pattern=$gap_string;}
												
												#-------------- now change identical bases with dots in $best_fitting_pattern
												$best_fitting_pattern=&change_bases_to_dots($best_fitting_pattern,$model_repeat_bases);
												
												
												my $gapless_best_fitting_pattern=$best_fitting_pattern;$gapless_best_fitting_pattern=~ s/-//g;
												my $no_of_bases_to_chop=length($gapless_best_fitting_pattern)+length($insertion_base);
												my $chopping_pattern=substr($right_flank,0,$no_of_bases_to_chop);
												
												#print "\$chopping_pattern=$chopping_pattern with \$no_of_bases_to_chop=$no_of_bases_to_chop\n";
												
												
												
												
												
												
												

												
												#--- now modify the previous repeat
												
												#$$modified_array[$#{$modified_array}]=~ s/$gapless_best_fitting_pattern$//;
												$right_flank=~ s/^$chopping_pattern//;
												
												#$r_start=$r_start-length($gapless_best_fitting_pattern);
												
												if($insertion_base ne "")
													{
														$position=$r_start+length($r_seq)-$no_of_trailing_gaps+$position;
														
														if($existing_insertion_bases)
															{
																$insertion_base=$insertion_base.",".$existing_insertion_bases;
																$position=$position.",".$existing_insertion_positions;
															}
															
														$comment="$insertion_base [$position]";
													}
												elsif($existing_insertion_bases)
													{
														
														$comment="$existing_insertion_bases [$existing_insertion_positions]";
													}	
													
												$r_seq=~ s/\-+$/$best_fitting_pattern/;
												
												#print "\nModified \$right_flank : $right_flank\n\n";
												
											}
										#elsif(length($s_seq1)>$avg_spacer_length*0.9)
										else{
											
											
												
												my $factor_1=$no_of_trailing_gaps+1;
												my $model_repeat_bases=substr($model_repeat,-$no_of_trailing_gaps); # this one will be one base short from the $bases to allow 1 base insertion
																								
											
												#-------- check for base(s) deletion and skip if deletion is true -------------										
										
												#if( (($r_length-$no_of_trailing_gaps)<($r_length*3/4)) and ($s_length<($avg_spacer_length*3/4)))
												if($s_seq=~/^$model_repeat_bases/)
													{
														#
													}
												elsif(($s_seq!~/[ATGC]/) or ( ($no_of_trailing_gaps>=($r_length*10)/100) and ($s_length<($avg_spacer_length-$avg_spacer_length*10/100))))
													{
														my $tmp_pos=$r_start+($r_length-$no_of_trailing_gaps);
														$comment="Deletion [$tmp_pos]";
														my $new_rec_line="$r_start\t$r_seq\t$s_seq\t$comment";
														
														#print "\$new_rec_line=$new_rec_line\n";
														push(@{$modified_array},$new_rec_line);
														
														next;
													}
												
												#------------------------------------------------------------------------------	
												my $no_of_tries=0;
												while($best_fitting_pattern=~/-$/)
													{													
													
														$no_of_tries++;
														if($no_of_tries>10){last;}
														
														my $bases=substr($s_seq,0,$factor_1); # this one will be one base longer that the $no_of_trailing_gaps 
														#print "\nTaken $no_of_trailing_gaps ( $bases ) from the spacer : $s_seq\n"; #exit;
														
														($best_fitting_pattern,$insertion_base,$position)=&get_best_fitting_string($range,$accession,$no_of_trailing_gaps,$bases,$model_repeat_bases,"RIGHT");												
														#print "$best_fitting_pattern,$insertion_base,$position\n\n";
														
														if($best_fitting_pattern eq ""){last;}
														
														my $no_of_t_gaps=0;
														if($best_fitting_pattern=~ /(\-+)$/)
															{
																#my $no_of_l_gaps=$1;
																
																#print "\$1=$1\n";
																$no_of_t_gaps =length($1);
															}
														#if($no_of_t_gaps>0 and length($insertion_base)>$no_of_t_gaps)
														#	{
														#		$best_fitting_pattern=substr($s_seq,0,$no_of_trailing_gaps);
														#		$insertion_base="";
														#		undef $position;
														#		#print "$best_fitting_pattern,$insertion_base,$position\n\n";
														#		
														#		last;
														#	}	
															
														#$factor_1=$factor_1+$no_of_t_gaps;
														if(length($s_seq)>$avg_spacer_length)
															{
																$factor_1=$factor_1+(length($s_seq)-$avg_spacer_length);
															}
														else{
																if($no_of_t_gaps)
																	{
																		$factor_1=$factor_1+$no_of_t_gaps;
																	}
																else{
																		$factor_1=$factor_1+1;
																	}		
															}
														if(($factor_1>(length($gap_string)*4)) or $factor_1>$s_length)
															{
																$best_fitting_pattern=substr($s_seq,0,$no_of_trailing_gaps);
																$insertion_base="";
																undef $position;
																last;
															}
														
													}
												
												#print "\$best_fitting_pattern=$best_fitting_pattern\n";
												
												if($best_fitting_pattern eq ""){$best_fitting_pattern=$gap_string;}
												
												#-------------- now change identical bases with dots in $best_fitting_pattern
												$best_fitting_pattern=&change_bases_to_dots($best_fitting_pattern,$model_repeat_bases);
												
												
												#---now get the total length of ungapped $gapless_best_fitting_pattern + the length of Insertion bases, grab the substr using this no. and chop off from the flank/spacer
												my $gapless_best_fitting_pattern=$best_fitting_pattern;$gapless_best_fitting_pattern=~ s/-//g;
												my $no_of_bases_to_chop=length($gapless_best_fitting_pattern)+length($insertion_base);
												my $chopping_pattern=substr($s_seq,0,$no_of_bases_to_chop);
												
												
												
												#--- now modify the spacer and repeat
												
												$s_seq=~ s/^$chopping_pattern//;
												
												#$r_start=$r_start-length($gapless_best_fitting_pattern);
												
												if($insertion_base ne "")
													{
														$position=$r_start+length($r_seq)-$no_of_trailing_gaps+$position;
														if($existing_insertion_bases)
															{
																$insertion_base=$insertion_base.",".$existing_insertion_bases;
																$position=$position.",".$existing_insertion_positions;
															}
															
														$comment="$insertion_base [$position]";
													}
												elsif($existing_insertion_bases)
													{
														
														$comment="$existing_insertion_bases [$existing_insertion_positions]";
													}	
													
												$r_seq=~ s/\-+$/$best_fitting_pattern/;
												
												#print "modified Repeat: $r_seq \tSpacer : $s_seq\n";
											}
										
										#---- now replace the leading gaps with bases from the flank/previous spacer
										
									}								

								$case_found=1;
									
							}
							
						#-----		
						
						my $repeat_seq;
						my $spacer_seq;
												
						#$repeat_seq=substr($species_seq,($r_start-1),$r_length);
						
						
						#---- now create two seq files and get the alignment
						my $region_len=$r_length+$s_length;
						my $r_stop=$r_start+$region_len;

						#print "$r_start\t$r_seq\t$s_seq\n";
						
						
						
						my $new_rec_line="$r_start\t$r_seq\t$s_seq\t$comment";
						push(@{$modified_array},$new_rec_line);
						#system("echo ");
						
			
					} 
		
		
		#-----update the @current array with modified array if case found
		#if($case_found==1)
		#	{
		#		push(@$modified_array,$$current_array[$#current_array]); 
		#		@{$current_array}=@{$modified_array};
		#	}
		
		return ($case_found);		
	}


sub get_best_fitting_string()
	{
		my($range,$accession,$no_of_gaps,$bases,$model_repeat_bases,$side)=@_;
		
		#print "\nInside sub: $no_of_gaps,$bases,$model_repeat_bases,$side\n\n"; #return;
		
		my $allowed_percent_similarity=33;
		
		my $best_fitting_pattern="";
		#my $model_repeat_bases="";
		my $insertion_base="";
		my $insertion_base_position="";
		
		my @array_bases;
		my @array_model_repeat_bases;
		
		$bases=~s/-//g;		
		if(not $bases){return ($best_fitting_pattern,$insertion_base,$insertion_base_position);}
		
		
		@array_bases=split('',$bases);	
		@array_model_repeat_bases=split('',$model_repeat_bases);
		
		#---- now get pairwise alignments using water
		if($no_of_gaps>1)
			{

						
				#for(my $i=$no_of_gaps;$i>0;$i--)
				#	{
					my $top_line="";
					my $bottom_line="--";
					
					my $n_o_g=0;
					my $difference_in_start_position=0;
					my $gap_increase=0;
					
					my $best_top_line="";
					my $best_bottom_line="";
					
					my $no_of_tries=0;
					my %hash_of_alignments;
					while($no_of_tries<3)
					{
						$no_of_tries++;
						
						my $gapopen=5+$gap_increase;
						my $gapextend=$gapopen-0.05; 						
						if($gapextend>10){$gapextend=9.5;}
						#system("echo 'Inside whhile2' >>log1.txt");											
						my $outfile=&run_water($range,$accession,$bases,$model_repeat_bases,$gapopen,$gapextend);					
						#---- now open the output file and get the alignment, store the alignment in a compiled file
						if(-e "$tmp_dir\/$outfile")
							{
								open(RD1,"$tmp_dir\/$outfile") or print "";#Error in module:get_best_fitting_string Line: __LINE__ \t$!";
								my @arr_rd=<RD1>;
								close(RD1);
								unlink("$tmp_dir\/$outfile");
								
								my $tl_count=0;
								my $bl_count=0;
								foreach my $line(@arr_rd)
									{
										#print $line,"\n";
										
										#system("echo '$line' >>log1.txt");
										
										chomp $line;$line=~s/\r//g;
										if(not $line or $line=~/#/){next;}
										elsif($line=~/\d+/)
											{
												#system("echo '$line' >>log1.txt");
												
												$line=~s/^\s+//;$line=~s/\s+/\t/g;
												my($start,$seq,$stop)=split('\t',$line);
												if($tl_count==0){$top_line=$line;$tl_count++;}
												elsif($bl_count==0){$bottom_line=$line;$bl_count++;}
											}	
										
										#if($line=~/score/i)
										#	{
										#		print "$line\n";
										#	}
									}
								
								#print "\n";
								
								
							}
						#else{next;}
						

						if($top_line eq "" or $bottom_line eq ""){next;}
						
						my($top_start,$top_seq,$top_stop)=split('\t',$top_line);
						my($bottom_start,$bottom_seq,$bottom_stop)=split('\t',$bottom_line);
						
						if(not defined $top_seq or not defined $bottom_seq){next;}
						
						$n_o_g=$bottom_line=~ s/-/-/g;
						if(not $n_o_g){$n_o_g=0;}
						$difference_in_start_position=$top_start-$bottom_start;
						
						
						#------ now calculate the score and store in an array
						my $similarity_score=&get_similarity_score($top_seq,$bottom_seq);
						$similarity_score=$similarity_score;#-$n_o_g;
						
						#print "\$similarity_score=$similarity_score\t\$n_o_g=$n_o_g\n\n";
						
						
						$hash_of_alignments{$similarity_score}{$n_o_g}=$top_line.",".$bottom_line;						
						
						if($n_o_g<=0 )      # max one insertion is allowed, ($top_start-$bottom_start)<=1 will achieve that
							{
								$best_top_line=$top_line;
								$best_bottom_line=$bottom_line;
								last;
							}
						elsif($no_of_tries>=3)
							{
								
								#---- now get the best scorring string---------
								my $best_rec_index=0;
								my $best_score=0;
								my $best_rec_line="";
								my $best_rec_nog=0;
								foreach my $s_score(sort{$b<=>$a}keys %hash_of_alignments)
									{
										foreach my $gaps(sort{$a<=>$b}keys %{$hash_of_alignments{$s_score}})
											{
												#print "$best_rec : $s_score\t$gaps\t$hash_of_alignments{$s_score}{$gaps}\n";
												
												
												$best_score=$s_score;
												$best_rec_nog=$gaps;
												$best_rec_line=$hash_of_alignments{$s_score}{$gaps};
												
												$best_rec_index++;
												last;
											}
										if($best_rec_index>0){last;}	
									}
								#----------------------------------------------
								#system("echo '$best_rec_line' >>log1.txt");
								if($best_rec_index>0) 
									{
										($best_top_line,$best_bottom_line)=split(',',$best_rec_line);
										
										if($best_bottom_line=~/-/){return("","","");}
										
									}
								else{	
										
										return("","","");
										#-------------------------------------------------------
										#$best_top_line="0\t$bases\t$no_of_gaps";
										#$best_bottom_line="0\t$model_repeat_bases\t$no_of_gaps";
										#		#$best_bottom_line="-";
										#		
										#$n_o_g=0;
										#$difference_in_start_position=0;
										#last;		
									}
								
								
								
								#$best_top_line="0\t$bases\t$no_of_gaps";
								#$best_bottom_line="0\t$model_repeat_bases\t$no_of_gaps";
								
								#$n_o_g=0;
								#$difference_in_start_position=0;
								#last;
								
							}	
						#$no_of_tries++;
						$gap_increase++;
					}	
						
						
					#print "\n\nBEST TOP LINE: $best_top_line\nBest Bottom line: $best_bottom_line\n";	
					
					

					
					#---- read this first
					#-- not working with ($top_start-$bottom_start)<=1 ; may need to individually delete each base at a time, and check the alignment
					#---- read above line
					#--------------------------------------------------
					
					
					#----- now complete check the alignment and add bases which is out side the alignment
					
					# Case 1: no gaps in the bottom line
					# case 2: one gap in the bottom line	
					# case 3: no alignment found	
					#system("echo 'TOP:$best_top_line' >>log1.txt");
					if(not defined $best_top_line or not defined $best_bottom_line){return("","","","","","");}
					
					my($top_start,$top_seq,$top_stop)=split('\t',$best_top_line);
					if(not defined $top_start or not defined $top_seq or not defined $top_stop){return("","","","","","");}
					$top_start=$top_start-1;
					$top_stop=$top_stop-1;
					#system("echo 'BOTTOM:$best_bottom_line' >>log1.txt");
					my($bottom_start,$bottom_seq,$bottom_stop)=split('\t',$best_bottom_line);
					if(not defined $bottom_start or not defined $bottom_seq or not defined $bottom_stop){return("","","","","","");}
					$bottom_start=$bottom_start-1;
					$bottom_stop=$bottom_stop-1;
					
					#if(not defined $top_seq or not defined $bottom_seq){return("","","","","","");}
					
					#print "\$difference_in_start_position=$difference_in_start_position\n";
					
					
					#--- read this first:
					#--- make the string in such way that there will be Gaps either end and finally just use the model string's length to use substr
					#-- for LEFT and RIGHT the conditions will vary
					# case 1: the matching bases can be exactly in the middle and there can be unmatched bases in both side
					# case 2: matching bases are on the left
					# case 3: matching bases are on the right
					# case 4: Insertion base(s) are the bases which have  higher position index than the length of the model repeat string after gap insertion, and if a natural insertion would have occured in water alignment
					
					my @arr_bfp=split("",$top_seq);	
					
					#---------- extend on the left -------------------------------------------									
					my $index_1=1;
					my @tmp_arr1;
					for(my $k1=$bottom_start-1;$k1>=0;$k1--)
						{
							if(($top_start-$index_1)>=0)
								{
									unshift(@tmp_arr1,$array_bases[$top_start-$index_1]);
								}
							else{
									unshift(@tmp_arr1,"-");
								}	
							
							$index_1++;
						}
									
					#---------- extend on the right -------------------------------------------				
					my $index_2=1;
					my @tmp_arr2;
									
					#print "\$bottom_stop=$bottom_stop \$#array_model_repeat_bases=$#array_model_repeat_bases\n";								
									
					for(my $k2=$bottom_stop+1;$k2<=$#array_model_repeat_bases;$k2++)
						{
							#print "Inside \$array_bases[$top_stop+$index_2]=$array_bases[$top_stop+$index_2]\n";
							if(($top_stop+$index_2)<=$#array_bases)
								{
									push(@tmp_arr2,$array_bases[$top_stop+$index_2]);
								}			
							else{
									push(@tmp_arr2,"-");
								}
							$index_2++;
						}
					if(defined $tmp_arr1[0])
						{						
							unshift(@arr_bfp,join("",@tmp_arr1));
						}
					if(defined $tmp_arr2[0])
						{		
							push(@arr_bfp,join("",@tmp_arr2));
						}	
									
										
					$best_fitting_pattern=join("",@arr_bfp);	
					#print "\$best_fitting_pattern=$best_fitting_pattern\n";	
					
					#system("echo 'BFS:$best_fitting_pattern' >>log1.txt");
					
					if($side eq "LEFT") #  all operations including base deletion should start on the left side
						{	
									#print join("",@array_bases),"\n";
									#---------- now get the insertion bases and position
									my $gap_less_bfp=$best_fitting_pattern;
									$gap_less_bfp=~s/-//g;
									
									
									if($bases=~s/$gap_less_bfp$//)
										{
											$insertion_base="";
											$insertion_base_position=0;
										}
									else{
											#print "inside\n";
											my @tmp_arr3=split($gap_less_bfp,$bases);
											
											$insertion_base=$tmp_arr3[1];
											$insertion_base_position=length($gap_less_bfp);
										}
									#$insertion_base=substr($bases,0,abs($#array_bases-$#array_model_repeat_bases));
									#$insertion_base_position=0;
									
									
									
									#print "\nA:	\$best_fitting_pattern=$best_fitting_pattern with $insertion_base\t found at $insertion_base_position\n";	
								
						
								
						}
					
					elsif($side eq "RIGHT") #  all operations including base deletion should start on the right side
						{				
									
									#print join("",@array_bases),"\n";
									#---------- now get the insertion bases and position
									my $gap_less_bfp=$best_fitting_pattern;
									$gap_less_bfp=~s/-//g;
									
									
									if($bases=~s/^$gap_less_bfp//)
										{
											$insertion_base="";
											$insertion_base_position=0;
										}
									else{
											#print "inside\n";
											my @tmp_arr3=split($gap_less_bfp,$bases);
											
											$insertion_base=$tmp_arr3[0];
											$insertion_base_position=0;
										}
									#$insertion_base=substr($bases,0,abs($#array_bases-$#array_model_repeat_bases));
									#$insertion_base_position=0;
									
									
									
									#print "\nB:	\$best_fitting_pattern=$best_fitting_pattern with $insertion_base\t found at $insertion_base_position\n";	
								
						}
					
			}
		else{
				if($side eq "LEFT")
					{	
						$best_fitting_pattern=$array_bases[$#array_bases];							
					}
				elsif($side eq "RIGHT")
					{
						$best_fitting_pattern=$array_bases[0];
					}
				
			}	
		
		return ($best_fitting_pattern,$insertion_base,$insertion_base_position);
	}







sub trim_repeat_ends_and_improve_score()
	{
		my($range,$minimum_repeat_length,$user_side_to_trim,$trimming_cutoff,$accession,$model_repeat,$avg_spacer_length,$current_array,$modified_array)=@_;
		#print "Going to trim the repeats for $accession:\n\n";
		
				
		
		my $case_found=0;
		my $left_flank="";
		my $right_flank="";
		
		#print "Inside:\t trim_repeat_ends_and_ improve_score\n";

		
		
		
		#---------- now open the sequence file and get the sequence string -------------
		open(SEQ,"$tmp_dir\/$accession\.fna") or print "$!";
		my @arr_seq=<SEQ>;
		close(SEQ);
		my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;
		
		
		
		
		#---- logic ------------------------------------------------------------------------
		#- Step 1: obtain the negative score of the array
		#- Step 2: determine the side(s) to trim
		#- Step 3: re-evaluate the score of the modified array, 
		#- Step 4: also check the spacer length distribution 
		
		
		#---------- obtain the overall array negative score and side specific degeneracy scores
		my @arr_repeats;
		my $no_of_repeats_with_no_degeneracy=0;
		
		my $array_score=0;
		my $left_degeneracy_score=0;
		my $center_degeneracy_score=0;
		my $right_degeneracy_score=0;
		
		
		#logic: devide the model repeat into 3 portions and check the degeneracy
		my $devision_length=int(length($model_repeat)/3);
		my $no_of_repeats=0;
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
						#print "@{$current_array-[0]}\n";
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
							
				#print "A: $current_line\n";	#next;			
				
				
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
						
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1];
				if(not defined $tmp_array[2])
					{
						$s_seq="";
					}
				else{
						$s_seq=$tmp_array[2];
					}		
																				
				$r_length=length($r_seq);
				#$s_length=length($s_seq);
				
				push(@arr_repeats,$r_seq);		
				
				if($r_seq!~/[ACGTU-]/i){$no_of_repeats_with_no_degeneracy++}
				
				#if($k1>4 and $k1<($#{$current_array}-1)) #----leaving the first and last repeats for future extension of array
				#	{
						my($lds,$cds,$rds)=&obtain_array_scores($r_seq,$model_repeat,$devision_length);		
						#print "$r_seq,$model_repeat,$devision_length\n\n";
						#$array_score=$array_score+$as;
						$left_degeneracy_score=$left_degeneracy_score+$lds;
						$center_degeneracy_score=$center_degeneracy_score+$cds;
						$right_degeneracy_score=$right_degeneracy_score+$rds;
				#	}
				
				$no_of_repeats++;
			}
		$array_score=$left_degeneracy_score+$center_degeneracy_score+$right_degeneracy_score;	
		#print "\n\$array_score=$array_score, $left_degeneracy_score,$center_degeneracy_score,$right_degeneracy_score\t\$no_of_repeats_with_no_degeneracy=$no_of_repeats_with_no_degeneracy\n";
			
		#------------- calculate the trimming_cutoff if set to AUTO or 0 ------------------------------
		if($trimming_cutoff=~/AUTO/ or $trimming_cutoff=~/^0/)
			{				
				#$trimming_cutoff=int(100/($no_of_repeats+1));	#;60-($no_of_repeats)*10;							
				
				#if($trimming_cutoff<20){$trimming_cutoff=20;} #--- setting the minimum boundary as for long arrays it will be negative value
				
				if($no_of_repeats<3)
					{
						$trimming_cutoff=50;						
					}
				elsif($no_of_repeats>=3 and $no_of_repeats <=6)
					{
						$trimming_cutoff=25;						
					}					
				else{
						$trimming_cutoff=33;						
					}
				
				
			}
		#----------------------------------------------------------------------------------------------
		
		
		my $trimming_factor=1;
		
		if($right_degeneracy_score<=$left_degeneracy_score)
			{
				$trimming_factor=$right_degeneracy_score;
			}
		elsif($left_degeneracy_score<$right_degeneracy_score)
			{
				$trimming_factor=$left_degeneracy_score;
			}	
		
		#---- check to see if it's worth to trim if $user_side_to_trim=~/NA/
		if($user_side_to_trim=~/NA/)
			{	
				if(length($model_repeat)<=$minimum_repeat_length or $array_score>=-1 or $trimming_factor==0 or $no_of_repeats_with_no_degeneracy >=int(($#arr_repeats+1)*(100-$trimming_cutoff)/100) )#or abs($trimming_factor)<$no_of_repeats_with_no_degeneracy)
					{
						#print "A.skipping with \$trimming_factor=$trimming_factor...\t";
						
						if( $no_of_repeats_with_no_degeneracy >=int(($#arr_repeats+1)*(100-$trimming_cutoff)/100))
							{
								#print "$no_of_repeats_with_no_degeneracy >=int(($#arr_repeats+1)*80/100) \n";
							}
						elsif(abs($trimming_factor)<$no_of_repeats_with_no_degeneracy)
							{
								#print "abs($trimming_factor)<$no_of_repeats_with_no_degeneracy\n";
							}
						elsif($array_score>=-1)
							{
								#print "$array_score>=-1\n";
							}
									
						return($model_repeat,$case_found);
					}	
				elsif( abs($trimming_factor)<($#arr_repeats+1)*$trimming_cutoff/100 )# or abs($trimming_factor)< int($#arr_repeats/2))     # side degeneracy rate >40% for trimming
					{
						my $score1=abs($trimming_factor)/$#arr_repeats;
						#print "B. Skipping trimming with \$trimming_factor=$trimming_factor and Score: $score1 ...\n\n";
						return($model_repeat,$case_found);	
					}
			
			}
		
		#print "Going to trim....\n";
		
		#--------- now determine the side(s) to trim
		my $side_to_trim;
		my @array_of_sides_to_trim;
		
		#if($right_degeneracy_score<=$left_degeneracy_score)
		#	{
		#		$side_to_trim="RIGHT";
		#	}			
		#else{
		#		$side_to_trim="LEFT";				
		#	}	
			
		my $final_trimmed_model_repeat=$model_repeat;		
		
		my $u_base_index;
		my $user_trimmed_model_repeat="";
		
		if($user_side_to_trim !~ /NA/)
			{
				undef @array_of_sides_to_trim;
							#push(@array_of_sides_to_trim,$user_side_to_trim);
				
				my @tmp_arr_1;
				if($user_side_to_trim=~/,/){@tmp_arr_1=split(',',$user_side_to_trim);}
				else{
						push(@tmp_arr_1,$user_side_to_trim);
					}
					
				my $to_trim_on_left=0;
				my $to_trim_on_right=0;
				
				foreach my $side_and_bases(@tmp_arr_1)
					{
						if($side_and_bases=~/RIGHT--/){next;}   # example RIGHT--1
						my $mr_length=0;
						if(defined $model_repeat and $model_repeat ne "")
							{
								$mr_length=length($model_repeat);
							}
								
						my($side,$bases)=split("-",$side_and_bases);
						#print "$mr_length - $to_trim_on_left - $bases : $user_side_to_trim\n";
						if($side =~/LEFT/ and int($mr_length - $to_trim_on_right - $bases)>=11)
							{
								$final_trimmed_model_repeat=substr($model_repeat,$bases);
								$u_base_index=$bases-1;
								
								push(@array_of_sides_to_trim,"LEFT");
								
								$to_trim_on_left=$bases;
							}
						elsif($side =~/RIGHT/ and int($mr_length - $to_trim_on_left - $bases)>=11)
							{
								$final_trimmed_model_repeat=substr($model_repeat,0,length($model_repeat)-$bases);
								$u_base_index=length($model_repeat)-$bases;
								
								push(@array_of_sides_to_trim,"RIGHT");
								
								$to_trim_on_right=$bases;
							}
					}		
			}	
		else{					
				if(abs($right_degeneracy_score)>0 and abs($right_degeneracy_score)>=($#arr_repeats+1)*$trimming_cutoff/100)
					{
						#$side_to_trim="RIGHT";
						push(@array_of_sides_to_trim,"RIGHT");
					}
				
				if(abs($left_degeneracy_score)>0 and abs($left_degeneracy_score)>=($#arr_repeats+1)*$trimming_cutoff/100)
					{
						#$side_to_trim="RIGHT";
						push(@array_of_sides_to_trim,"LEFT");
					}
			}
	

		#print "\@array_of_sides_to_trim=@array_of_sides_to_trim<br>\n\n";
		
		
		
		
		
		
		
		
		foreach my $side_to_trim(@array_of_sides_to_trim)
			{
				#print "\$side_to_trim=$side_to_trim\n";
				undef @arr_repeats;
			
				for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
					{
								#print "@{$current_array-[0]}\n";
						my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
									
								#print "A: $current_line\n";
									
						my @tmp_array=split('\t',$current_line);
						my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
								
						$r_start=$tmp_array[0];
						$r_seq=$tmp_array[1];
						$s_seq=$tmp_array[2];																		
						
						push(@arr_repeats,$r_seq);		
						
					}
				my($trimmed_model_repeat,$base_index);
				
				if($user_side_to_trim =~ /NA/)
					{
						($trimmed_model_repeat,$base_index)=&obtain_trimming_model_repeat_and_base_index($minimum_repeat_length,$trimming_cutoff,$side_to_trim,$devision_length,$final_trimmed_model_repeat,\@arr_repeats);
						#print "A: \$trimmed_model_repeat=$trimmed_model_repeat\t\$base_index=$base_index\n";
						#system("echo '\$trimmed_model_repeat=$trimmed_model_repeat\t\$base_index=$base_index' >log.log");
						if($base_index<0){next;}
					}
				else{
						$trimmed_model_repeat=$final_trimmed_model_repeat;
						$base_index=$u_base_index;
					}	
					
					
				$final_trimmed_model_repeat=$trimmed_model_repeat;
				
				#print "\$side_to_trim=$side_to_trim\t\$trimmed_model_repeat=$trimmed_model_repeat\t\$base_index=$base_index\n";	
				
				
				
				
				
				
				
				
				#print "\$side_to_trim=$side_to_trim with \$base_index=$base_index\n\n";
				
				#----------- now create a new array with the trimmed_model_repeat
				my @arr_model_repeat=split('',$model_repeat);
				
				for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
					{
								#print "@{$current_array-[0]}\n";
						my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
									
						#print "B: $current_line\n";
									
						my @tmp_array=split('\t',$current_line);
						my($r_start,$r_length,$s_length,$r_seq,$s_seq);
						my $comment;
						
						my $original_s_seq_stop;
						my $no_of_insertions=0;	
						
							
						$r_start=$tmp_array[0];
						$r_seq=$tmp_array[1];
						
						if(not defined $tmp_array[2])
							{
								$s_seq=""; 
							}
						else{
								$s_seq=$tmp_array[2]; 
							}
							
						my $gapless_s_seq=$s_seq; if(defined $gapless_s_seq){chomp $gapless_s_seq;$gapless_s_seq=~s/\-+//g;}
						
												
						my %hash_of_insertion_positions;
						my $comment_del_part="";
						if($tmp_array[3] and $tmp_array[3]!~/^Del/)
							{
								$comment=$tmp_array[3]; chomp $comment; $comment=~s/^\s+//;
								#if($comment=~/^Del/){next;}
								
								my @tmp_arr1=split(' ',$comment);
														
														#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
								my $insertion_bases=$tmp_arr1[0];
								my $insertion_positions=$tmp_arr1[1];
								#print "\t\$insertion_bases=$insertion_bases and \$insertion_positions=$insertion_positions\n";
								if($insertion_positions and $insertion_positions!~/Del/)
									{
									
									   $insertion_positions=~s/\[//g;
									   $insertion_positions=~s/\]//g;
									   #print "\$cur_insertion_bases=$cur_insertion_bases\n";
										
										my @tmp_arr2; my @tmp_arr3;
										
										if($insertion_bases=~/,/)
											{ 			
												@tmp_arr2=split(',',$insertion_bases);
											}
										else{
												push(@tmp_arr2,$insertion_bases);
											}	
											
										if($insertion_positions=~/,/)
											{
												@tmp_arr3=split(',',$insertion_positions);
											}
										else{			
												push(@tmp_arr3,$insertion_positions);
											}	
											
										for( my $p=0;$p<=$#tmp_arr3;$p++)
											{
												my $pos=$tmp_arr3[$p];
												if(defined $tmp_arr2[$p])
													{
														chomp $tmp_arr2[$p];$tmp_arr2[$p]=~s/^\s+//;
														if($tmp_arr2[$p] ne "")
															{
																$hash_of_insertion_positions{$pos}=$tmp_arr2[$p];
															}
													}
												
												
												#print "\tA:\$hash_of_insertion_positions{$pos}=$hash_of_insertion_positions{$pos}\n";
											}
									
									
																			
										$insertion_bases=~s/,//g;
										$no_of_insertions=length($insertion_bases);		
										
									
									}
							
								#--------------- save the Del part saved
								
								if(defined $tmp_arr1[2] and $tmp_arr1[2]=~/Del/)
									{
										$comment_del_part=$tmp_arr1[2]." ".$tmp_arr1[3];
									}
							}
						else{
								if(defined $comment and $comment=~/^Del/)
									{
										$comment_del_part=$comment;
									}
								
								$comment="";
							}		
						
						my $gapless_r_seq=$r_seq;$gapless_r_seq=~s/-//g;
						$original_s_seq_stop=$r_start+length($gapless_r_seq)+$no_of_insertions+length($gapless_s_seq);		
																				
						#print "\$original_s_seq_stop=$original_s_seq_stop\t\t";
						
						#----- get the previous record ---------------------------------------------------
						my $p_no_of_insertions=0;
						my($p_r_start,$p_r_seq,$p_s_seq,$p_comment);#=split('\t',$$current_array[$k1-1]);
						
						my @tmp_array2=split('\t',$$current_array[$k1-1]);
						$p_r_start=$tmp_array2[0]; 
						$p_r_seq=$tmp_array2[1];
						$p_s_seq=$tmp_array2[2];
						
						if($tmp_array2[3])
							{
								$p_comment=$tmp_array2[3];chomp $p_comment;$p_comment=~s/^\s+//;
								
								if($p_comment ne "" and $p_comment!~/^Del/)
									{
										#$comment=$tmp_array[3];
										
										#------ if there are any insertion bases whose position is out side the current repeat, then put them accordingly in the spacer 
										my @tmp_arr1=split(' ',$p_comment);
														
														#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
										my $p_insertion_bases=$tmp_arr1[0];
										my $p_insertion_positions=$tmp_arr1[1];
										
										if($p_insertion_positions)
										{
										
										$p_insertion_positions=~s/\[//g;
										$p_insertion_positions=~s/\]//g;
																#print "\$cur_insertion_bases=$cur_insertion_bases\n";
												
																		
										$p_insertion_bases=~s/,//g;
										$p_no_of_insertions=length($p_insertion_bases);		
												
										#$p_comment=$comp_p_insertion_bases." [".$p_insertion_positions."]"; 
										}
									}
							}
						else{
								$p_comment="";
							}
						
						
						
						
						my $new_repeat="";
						my $new_spacer="";
						my $trimmed_repeat="-";
						
						my @arr_current_repeat=split('',$r_seq);
						
						
						
						
						if($side_to_trim eq "RIGHT")
							{
								$new_repeat=substr($r_seq,0,$base_index);
								
								
								
								#----- check if any Insertion(s) was there before which falls outside
								my $gapless_new_repeat=$new_repeat;$gapless_new_repeat=~s/-//g;
								foreach my $pos(sort{$a<=>$b} keys %hash_of_insertion_positions)
									{
										if(not defined $hash_of_insertion_positions{$pos}){next;}
										
										my $t_r_stop=$r_start+length($gapless_new_repeat);
										if($pos>=($r_start + length($gapless_new_repeat) - (length($hash_of_insertion_positions{$pos})-1) ))
											{
												delete $hash_of_insertion_positions{$pos};
											}
										else{	
												#print "\$t_r_stop=$t_r_stop\t\$hash_of_insertion_positions{$pos}=$hash_of_insertion_positions{$pos}\n";	
											}	
									}										
							}
							
								#------------ for LEFT -----------------------------------------------
						else{ 
								#-------------- get the $no_of_leading_gaps first ----------------------
								my $no_of_leading_gaps=0;
								if($r_seq=~ /^(\-+)/)
									{
										$no_of_leading_gaps =length($1);								
									}
									
								if($no_of_leading_gaps>0)
									{
										#print "\t$r_start\n";
										if($no_of_leading_gaps>$base_index+1)
											{
												$r_start=$r_start;
											}
										elsif($base_index+1>=$no_of_leading_gaps)
											{
												$r_start=$r_start+$base_index+1-$no_of_leading_gaps;
											}
										#print "\t$r_start\n";		
									}
								else{	
										#------ if the gaps are in the middle ---------------------------------
										my $no_of_total_gaps_in_chopping_string=0;
										my @arr_r_seq=split('',$r_seq);
										for(my $z=0;$z<=$base_index;$z++)
											{
												if($arr_r_seq[$z] eq "-"){$no_of_total_gaps_in_chopping_string++;}
											}
											
										$r_start=$r_start+$base_index+1-$no_of_total_gaps_in_chopping_string;			#length($trimmed_repeat);	
									}									
										
							
								$new_repeat=substr($r_seq,($base_index+1));
								
								#----- check if any Insertion(s) was there before which falls outside
								#my $gapless_new_repeat=$new_repeat;$gapless_new_repeat=~s/-//g;
								
								#---- the following block will progressively increase the $r_start
								foreach my $pos(sort{$a<=>$b} keys %hash_of_insertion_positions)
									{
										#if($pos<=($r_start+($base_index+1)))      #---- as $r_start still pointing to the original $r_start
										if($pos<=($r_start+0))      #---- as $r_start now pointing to the corrected start
											{
												if(not defined $hash_of_insertion_positions{$pos}){next;}
												
												$r_start=$r_start+length($hash_of_insertion_positions{$pos}); #--- as there could be multiple bases per insertion
												
												#print "B:\$hash_of_insertion_positions{$pos}=$hash_of_insertion_positions{$pos}\n";
												delete $hash_of_insertion_positions{$pos};
											}
									}
								

								
								if($k1>40000)  #-- not required, deleate after checking ---
									{	
										#print "\tA:$p_r_start\n";
										#$p_r_start=$p_r_start+$base_index+1;
										#print "\tB:$p_r_start\n";							
										#$p_s_seq=$p_s_seq.$trimmed_repeat;
										#
										#---- here get the $p_s_seq from species_seq 
										my $gapless_p_r_seq=$p_r_seq; $gapless_p_r_seq=~s/-//g;
										my $gapless_p_s_seq=$p_s_seq; $gapless_p_s_seq=~s/-//g;
										
										my $p_s_seq_start;
										
										if($k1==5)  #--- this is a raw fix, working but check later
											{
												$p_s_seq_start=$p_r_start+length($gapless_p_r_seq)-1 + $p_no_of_insertions + $base_index+1;
											}
										else{
												$p_s_seq_start=$p_r_start+length($gapless_p_r_seq)-1 + $p_no_of_insertions;
											}	
										
										
										
										my $p_s_seq_stop=$r_start-1;
										
										$p_s_seq=substr($species_seq,$p_s_seq_start,length($gapless_p_s_seq)+ $base_index+1);	
										
										if(not $p_s_seq or $p_s_seq eq ""){$p_s_seq="-";}
										#print "$k1:$p_r_start\t$p_r_seq\t$p_s_seq\t$p_comment\n";
										
										$$current_array[$k1-1]="$p_r_start\t$p_r_seq\t$p_s_seq\t$p_comment";
									}		

								 	
							}	
						
						
						
						
						
						
						
						
						#----- now recreate the $comment
						my $remaining_no_of_insertions=keys %hash_of_insertion_positions;
						my $gapless_new_repeat=$new_repeat;$gapless_new_repeat=~s/-//g;
						if($remaining_no_of_insertions>0)
							{
								my $c_insertion_bases="";
								my $c_insertion_positions="";
								
								
								foreach my $pos(sort{$a<=>$b} keys %hash_of_insertion_positions)
									{
										#print "\$hash_of_insertion_positions{$pos}=$hash_of_insertion_positions{$pos}\n";
										
										if($pos>$r_start and $pos<($r_start+length($gapless_new_repeat)+$remaining_no_of_insertions))
											{
												$c_insertion_bases=$c_insertion_bases.",".$hash_of_insertion_positions{$pos};
												$c_insertion_positions=$c_insertion_positions.",".$pos;
												
											}	
										else{
												delete $hash_of_insertion_positions{$pos};
											}	
									}
								
								$c_insertion_bases		=~s/^,//;										
								$c_insertion_positions	=~s/^,//;
								
									
									
								#if($tmp_arr1[2] and $tmp_arr1[3]) #--- if deletion was present
								#	{				
								#		$comment=$c_insertion_bases." [".$c_insertion_positions."] $tmp_arr1[2] $tmp_arr1[3]"; 
								#	}
								#else{
								#		$comment=$c_insertion_bases." [".$c_insertion_positions."]"; 
								#	}		
								
								if($c_insertion_bases ne "")
									{
										#print "\$comment=$comment\n";
										$comment=$c_insertion_bases." [".$c_insertion_positions."]"; 
									}
							}
						else{
								if($comment_del_part)
									{
										$comment=$comment_del_part;
									}
								else{	
										$comment="";
									}	
							}
						
						
						
						
						
						#------------------------------------ get the spacer seq ------------------------------------------------
						my $remaining_insertion_bases=0;
						
						foreach my $pos(keys %hash_of_insertion_positions)
							{
								$remaining_insertion_bases=$remaining_insertion_bases+length($hash_of_insertion_positions{$pos});
							}
						
						my $c_s_seq_start=$r_start+length($gapless_new_repeat)+$remaining_insertion_bases -1;
						
						#print "\$c_s_seq_start=$c_s_seq_start\n";
						
						if($side_to_trim eq "LEFT")
							{
								#print "$k1:$$current_array[$k1+1]\n";
								#------- check if there are any insertions falling withing the trimming bases of next repeat, casue they need to be added to the current spacer (by checking the length of the insertions)
								my $no_of_ins_in_next_repeat=0;
								
								if($k1<$#{$current_array}-1)
									{
										my $next_line=$$current_array[$k1+1]; chomp $next_line; $next_line=~s/\r+//g; $next_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
										my @arr_nxt1=split('\t',$next_line);
										if(defined $arr_nxt1[3] and $arr_nxt1[3]!~/^Del/)
											{
												#print "$arr_nxt1[3]\n";
												my %tmp_hash_of_insertions;
												my $tmp_no_of_ins=&get_number_of_insertions($arr_nxt1[3],\%tmp_hash_of_insertions);
												
												#---- check number of insertions falling within the chopping bases
												my $nxt_repeat_start=$arr_nxt1[0];
												foreach my $pos(keys %tmp_hash_of_insertions)
													{
														if($pos<=($nxt_repeat_start+$base_index+1))
															{
																$no_of_ins_in_next_repeat=$no_of_ins_in_next_repeat+length($tmp_hash_of_insertions{$pos});
															}
													}
												
											}
									}
								
								if($c_s_seq_start<length($species_seq))
									{
										$new_spacer=substr($species_seq,$c_s_seq_start,($original_s_seq_stop-$c_s_seq_start-1)+$base_index+1 + $no_of_ins_in_next_repeat);
									}	
							}
						else{
								#---- for RIGHT
								
								
								
								
								
								
								if(length($species_seq)>=($c_s_seq_start+($original_s_seq_stop-$c_s_seq_start-1)))
									{
										$new_spacer=substr($species_seq,$c_s_seq_start,($original_s_seq_stop-$c_s_seq_start-1));
									}
								elsif(length($species_seq)>$c_s_seq_start)
									{
										$new_spacer=substr($species_seq,$c_s_seq_start);
									}
								else{
										$new_spacer="-";
									}			
							}	
						
						if($k1==$#{$current_array}-1){$new_spacer=$s_seq;}
						#---- now remove leading gaps from the spacer sequence 						
						if(not defined $new_spacer or $new_spacer eq ""){$new_spacer="-";}
						$new_spacer=~s/-+/-/g;
						#------ now update the record with the new repeats and spacers	
						$$current_array[$k1]="$r_start\t$new_repeat\t$new_spacer\t$comment";
						#print "\tC: $$current_array[$k1]\n";	
						
						
					}	# end of for $k1
				
				#----- check that user selected bases are trimmed already, and if so return------------------------------------
				
				#if($user_side_to_trim ne "NA" and length($model_repeat)>=(length($final_trimmed_model_repeat)+$user_no_of_bases_to_trim) )
				#	{
				#		#$model_repeat=$final_trimmed_model_repeat;	
				#		last;
				#	}
					
				#---end of while loop	
						
					
			}	
			
		if($model_repeat ne $final_trimmed_model_repeat){$case_found=1;}
		$model_repeat=$final_trimmed_model_repeat;	
			
			#print "\$model_repeat=$model_repeat\n\n";
				
			return($model_repeat,$case_found);		
	}


					
sub obtain_trimming_model_repeat_and_base_index()
	{
		my($minimum_repeat_length,$trimming_cutoff,$side_to_trim,$devision_length,$model_repeat,$arr_repeats)=@_;
		
		
		
		#my $region_degeneracy_cutoff=int($trimming_cutoff/2);
		
		#print "\n$minimum_repeat_length,$trimming_cutoff,$side_to_trim,$devision_length,$model_repeat,$arr_repeats\n\n";
		my $trimmed_model_repeat=$model_repeat;
		my $base_index=-1;
		
		my %hash_of_bases;
		my @arr_sorted_bases;
		#foreach my $repeat(@{$arr_repeats})
		#	{
		#		#print "$repeat\n";
		#		
		#		my @arr_1=split('',$repeat);
		
		
		#print "\n\nNo. of repeats:",$#{$arr_repeats},"\n";		
				
		if($side_to_trim eq "RIGHT") #----- for RIGHT side ----------------------------
			{
								
				foreach my $repeat(@{$arr_repeats})
					{
						#print "$repeat\n";
						
						my @arr_1=split('',$repeat);
						
						for(my $i=$#arr_1;$i>=(2*$devision_length);$i--)
							{
								if($arr_1[$i] ne ".")
									{
										if(not $hash_of_bases{$i})
											{
												$hash_of_bases{$i}=1;
											}
										else{
												$hash_of_bases{$i}=$hash_of_bases{$i}+1;
											}
									}
							}
					}	
						
				#my $highest_value=0;	b
				foreach my $baseIndex(sort{$a<=>$b}keys %hash_of_bases)
					{
						
						#if($hash_of_bases{$baseIndex}<$highest_value or $baseIndex<22){next;}
						my $base_occurrence=$hash_of_bases{$baseIndex};
						my $column_degeneracy=int(($base_occurrence/($#{$arr_repeats}+1))*100 );
						
						#if($hash_of_bases{$baseIndex}<=1 or $hash_of_bases{$baseIndex}<=int($#{$arr_repeats}*$trimming_cutoff/100))#or $hash_of_bases{$baseIndex} <= int($#{$arr_repeats}*$trimming_cutoff/100))
						
						
						
						if($base_occurrence<1 or $column_degeneracy < $trimming_cutoff or $baseIndex<11) #----------- hard stop at 11
							{
								#print "$baseIndex: $hash_of_bases{$baseIndex} $column_degeneracy <= $trimming_cutoff\n";
								next;
							}
						else{
								#print "$baseIndex: $hash_of_bases{$baseIndex} $column_degeneracy <= $trimming_cutoff\n";
							}
								
						if(defined $arr_sorted_bases[0])
							{
								#print "\$baseIndex=$baseIndex\n";
								if($baseIndex < $arr_sorted_bases[0])
								#if($baseIndex<$arr_sorted_bases[0] and $hash_of_bases{$baseIndex} > int($#{$arr_repeats}*0.5))								
									{
										unshift(@arr_sorted_bases,$baseIndex);	
									}
								else{
										push(@arr_sorted_bases,$baseIndex);
									}		
							}
						else{
								#if($hash_of_bases{$baseIndex} > int($#{$arr_repeats}*$trimming_cutoff/100))
								#	{
										push(@arr_sorted_bases,$baseIndex);
										#$highest_value=$hash_of_bases{$baseIndex};	
								#	}
									
							}		
					}

							
			}
		#------------ now for the LEFT side -------------	
		else{
				foreach my $repeat(@{$arr_repeats})
					{
						#print "$repeat\n";
						
						my @arr_1=split('',$repeat);
			
						for(my $i=0;$i<=$devision_length;$i++)
							{
								if($arr_1[$i] ne ".")
									{
										if(not $hash_of_bases{$i})
											{
												$hash_of_bases{$i}=1;
											}
										else{
												$hash_of_bases{$i}=$hash_of_bases{$i}+1;
											}
									}
							}
					}		
						
				my $highest_value;	
				foreach my $baseIndex(sort{$b<=>$a}keys %hash_of_bases)
					{
						
						my $base_occurrence=$hash_of_bases{$baseIndex};
						my $column_degeneracy=int( ($base_occurrence/($#{$arr_repeats}+1))*100 );
						
						#print "\$baseIndex=$baseIndex \$hash_of_bases{$baseIndex}=$hash_of_bases{$baseIndex}\n";
						
						#if((length($model_repeat)-$baseIndex)<$minimum_repeat_length or $hash_of_bases{$baseIndex}<=1 or $hash_of_bases{$baseIndex} <= int($#{$arr_repeats}*$trimming_cutoff/100))
						if($base_occurrence<1 or $column_degeneracy < $trimming_cutoff or (length($model_repeat)-$baseIndex < 11)) #----- hard stop at 11
							{
								#print "\$baseIndex=$baseIndex $base_identity >= $trimming_cutoff\n";
								#print "skiping with (length($model_repeat)-$baseIndex)<$minimum_repeat_length or $hash_of_bases{$baseIndex}<=1 or $hash_of_bases{$baseIndex} <= int($#{$arr_repeats}*$trimming_cutoff/100)........\n";
								next;
							}
						
						
						if(defined $arr_sorted_bases[0])
							{
								#print "\$baseIndex=$baseIndex\n";
								if($baseIndex>$arr_sorted_bases[0])
								#if($baseIndex<$arr_sorted_bases[0] and $hash_of_bases{$baseIndex} > int($#{$arr_repeats}*0.5))								
									{
										unshift(@arr_sorted_bases,$baseIndex);	
									}
								else{
										push(@arr_sorted_bases,$baseIndex);
									}		
							}
						else{
								#if($hash_of_bases{$baseIndex} > int($#{$arr_repeats}*$trimming_cutoff/100))
								#	{
										push(@arr_sorted_bases,$baseIndex);
										#$highest_value=$hash_of_bases{$baseIndex};	
								#	}
									
							}										
					}
						
					
			}	
			
		
		
		
		
		#print "\n\@arr_sorted_bases=@arr_sorted_bases for \$side_to_trim=$side_to_trim\n";
		
		#----------- now for each of these stored bases, calculate the region degeneracy score,
		if(not defined $arr_sorted_bases[0])
			{
				if($side_to_trim eq "LEFT")
					{
						$trimmed_model_repeat=$model_repeat;
						$base_index=-1;
					}
				else{
						$trimmed_model_repeat=$model_repeat;
						$base_index=length($model_repeat);
							
						#print "\$base_index=$base_index\n";
					}	
			}
		else{	
				#print "\n\@arr_sorted_bases=@arr_sorted_bases for \$side_to_trim=$side_to_trim\n";
				
				foreach my $b_index(@arr_sorted_bases)
					{
						#print "\n";
						my $number_of_degenerated_bases=0;
						my $total_number_of_bases=0;
						my $no_of_gaps_found=0;
						
						foreach my $c_repeat(@{$arr_repeats})
							{
								my $trimmed_part="";
								
								if($side_to_trim eq "LEFT")
									{
										$trimmed_part=substr($c_repeat,0,$b_index+1);
									}
								else{
										$trimmed_part=substr($c_repeat,$b_index);
										
									}	
										
								#		$total_number_of_bases=$total_number_of_bases+length($trimmed_part);
										
								#		$trimmed_part=~s/\.//g;
										
								#		$number_of_degenerated_bases=$number_of_degenerated_bases+length($trimmed_part);
										
										#print "\$trimmed_part $trimmed_part\n";
										
										
								$total_number_of_bases=$total_number_of_bases+length($trimmed_part);
										
								if($trimmed_part=~ /-/)
									{
										my $n_o_g= $trimmed_part=~s/-/-/g;
										   $no_of_gaps_found=$no_of_gaps_found+$n_o_g;
									}	
										
								$trimmed_part=~s/\.//g;
										
								$number_of_degenerated_bases=$number_of_degenerated_bases+length($trimmed_part);
										
										#print "\$trimmed_part=$trimmed_part\n";
										
							}
							
						my $region_degeneracy_score=($number_of_degenerated_bases/$total_number_of_bases)*100 + 2*$no_of_gaps_found; # --- or just *$trimming_cutoff
						
						#print "\$b_index=$b_index $region_degeneracy_score =($number_of_degenerated_bases/$total_number_of_bases)*100*($trimming_cutoff/100)\n with \$model_repeat=$model_repeat\n";
						
						#my $region_degeneracy_score=$number_of_degenerated_bases/($#{$arr_repeats}+1);
						#my $region_degeneracy_score=($number_of_degenerated_bases/$total_number_of_bases)*100; #int($trimming_cutoff/100);
						$region_degeneracy_score=sprintf("%.0f",$region_degeneracy_score);
						
						my $length_of_region=0;
						
						if($side_to_trim =~/LEFT/)
							{
								$length_of_region=$b_index+1;
							}
						else{
								$length_of_region=length($model_repeat)-$b_index+1;
							}	
						
						my $region_degeneracy_cutoff=$trimming_cutoff; #int($trimming_cutoff+$length_of_region);#sprintf("%.0f",($total_number_of_bases*$trimming_cutoff)/100);
						
						#print "$side_to_trim \$b_index=$b_index \$region_degeneracy_score=$region_degeneracy_score  \$region_degeneracy_cutoff=$region_degeneracy_cutoff  \$trimming_cutoff=$trimming_cutoff\n";
						
						if($region_degeneracy_score > $region_degeneracy_cutoff)
							{
								#print "\$b_index=$b_index: \$region_degeneracy_score $region_degeneracy_score\n";
								if($side_to_trim eq "LEFT")
									{
										
										$trimmed_model_repeat=substr($model_repeat,$b_index+1);	
										$base_index=$b_index;	
										last;
									}
								else{
										$trimmed_model_repeat=substr($model_repeat,0,$b_index);	
										$base_index=$b_index;
										last;
									}		
							}
							
						#print "Region degeneracy cutoff for $b_index= $region_degeneracy_cutoff with \$total_number_of_bases=$total_number_of_bases and \$number_of_degenerated_bases=$number_of_degenerated_bases\n";
							
					}
				
				
				
			}	
		
		
		#print "\$base_index=$base_index\n";
		#$base_index=$arr_sorted_bases[0];
		
		#print "$model_repeat\n$trimmed_model_repeat\n\n";	
		return($trimmed_model_repeat,$base_index);
	}





sub increase_all_repeat_lengths()
	{
		my($range,$user_side_to_increase_length,$user_no_of_bases_to_increase_length,$accession,$model_repeat,$current_array,$modified_array)=@_;
		
		my $case_found=0;
		if($user_side_to_increase_length=~/NA/)
			{				
				return($model_repeat,$case_found);
			}
		#print "Going to increase the repeat_lengths for an the array of $accession: with $user_side_to_increase_length,$user_no_of_bases_to_increase_length,$accession,$model_repeat\n\n";
		
		#---------- now open the sequence file and get the sequence string -------------
		open(SEQ,"$tmp_dir\/$accession\.fna") or print "$!";
		my @arr_seq=<SEQ>;
		close(SEQ);
		my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;
		#-------------------------------------------------------------------------------
		
		my %hash_of_repeats;
		my @arr_repeats;			
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
				#print "@{$current_array-[0]}\n";
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
							
				#print "C: $current_line\n";
							
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq);
				my $comment="";
				my $current_spacer_bases_to_chop="";
						
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1]; my $r_seq_1=$r_seq;$r_seq_1=~s/-//g;
				if(not defined $tmp_array[2]){$tmp_array[2]="";}
				$s_seq=$tmp_array[2]; my $s_seq_1=$s_seq;$s_seq_1=~s/-//g;
				
				my $no_of_insertions=0;
				
				if($tmp_array[3])
					{
						$comment=$tmp_array[3];	$comment=~s/^\s+//;
						if($comment)
							{					
								if($comment=~/^Del/)    # --- deletion anywhere ----
									{
										#$array_total_degeneracy=$array_total_degeneracy+1;
										#next;
									}
								else{										
										my @tmp_arr1=split(' ',$comment);
												
												#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
										my $insertion_bases=$tmp_arr1[0];
										#my $insertion_positions=$tmp_arr1[1];
										#	$insertion_positions=~s/\[//g;
										#	$insertion_positions=~s/\]//g;
														#print "\$cur_insertion_bases=$cur_insertion_bases\n";
										
										#my $comp_p_insertion_bases=$insertion_bases;
										#$comp_p_insertion_bases=~tr/ACGT/TGCA/;   #--- no need to reverse the string
										
																						
										$insertion_bases=~s/,//g;
										$no_of_insertions=length($insertion_bases);		
										
										#$comment=$comp_p_insertion_bases." [".$insertion_positions."]";
									}
							}
						
					}					
				
				
				#---------------------------------------------------------------------
				my ($p_r_start,$p_r_seq,$p_s_seq);
				my $p_comment="";
				
				my $previous_spacer_bases_to_chop="";
				

										
				if($k1<$#{$current_array}-1)
					{	
						if(length($s_seq_1)>=$user_no_of_bases_to_increase_length)
							{
								$current_spacer_bases_to_chop=substr($s_seq_1,0,$user_no_of_bases_to_increase_length);	
							}
						else{
								#--- add the available bases to $current_spacer_bases_to_chop, and then add - to make it equal length
								$current_spacer_bases_to_chop=$s_seq_1;
								for(my $i=1;$i<=$user_no_of_bases_to_increase_length-length($s_seq_1);$i++) #Remember, starting from 1
									{
										$current_spacer_bases_to_chop=$current_spacer_bases_to_chop."-";
									}
							}		
					}
				else{
						if((($r_start-1+length($r_seq_1)+$no_of_insertions)+$user_no_of_bases_to_increase_length)<=length($species_seq))
							{
								if(($r_start+length($r_seq_1)-1+$no_of_insertions)<length($species_seq))
									{
										$current_spacer_bases_to_chop=substr($species_seq,($r_start+length($r_seq_1)-1+$no_of_insertions),$user_no_of_bases_to_increase_length);
									}	
							}
						else{
								#----- take from 0 position
								#$current_spacer_bases_to_chop=substr($species_seq,0,$user_no_of_bases_to_increase_length);
								for(my $i=1;$i<=$user_no_of_bases_to_increase_length;$i++)
									{
										$current_spacer_bases_to_chop=$current_spacer_bases_to_chop."-";
									}
							}	
					}		


						
				#-------------------------------------------------------------------
				#my $tmp_seq1=substr($species_seq,2664095,210);
				#print "\$tmp_seq1=$tmp_seq1\n\n";
				
				my $new_r_seq;
				
				#my $next_spacer_seq="";
				#my $next_spacer_bases_to_chop="";
				
				if($user_side_to_increase_length eq "RIGHT" and $user_no_of_bases_to_increase_length>0)
					{
						#----- if the repeat ends with gaps, then no bases should be taken from spacer
						if($r_seq=~/-$/)
							{
								$current_spacer_bases_to_chop="";
								for(my $i=0;$i<$user_no_of_bases_to_increase_length;$i++)
									{
										$current_spacer_bases_to_chop=$current_spacer_bases_to_chop."-";
									}
							}
						
						my $modified_r_start=$r_start;
						
						my $old_r_seq=&change_dots_to_bases($r_seq,$model_repeat);
						
						my $modified_repeat=$old_r_seq.$current_spacer_bases_to_chop;
						my $modified_spacer=$s_seq_1;
						
						#--- remove the "-"s from $current_spacer_bases_to_chop
						$current_spacer_bases_to_chop=~s/-//g;
						$modified_spacer=~s/^$current_spacer_bases_to_chop//;
						
						if($modified_spacer eq ""){$modified_spacer="-";}
						$s_seq=$modified_spacer;$s_seq_1=$s_seq;$s_seq_1=~s/-//g;
						
						my $modified_rec_line="";
						
						if($k1<$#{$current_array}-1)
							{
								$modified_rec_line="$r_start\t$modified_repeat\t$modified_spacer\t$comment";
							}
						else{
								$modified_rec_line="$r_start\t$modified_repeat\t|\t$comment";
							}		
						#print "\t$modified_rec_line\n";
						push(@{$modified_array},$modified_rec_line);
						
						$new_r_seq=$modified_repeat;
					}
				
				elsif($user_side_to_increase_length eq "LEFT" and $user_no_of_bases_to_increase_length>0)
					{					
						
						if($k1==4)
							{
								if(($r_start-1-$user_no_of_bases_to_increase_length)>=0 and ($r_start-1-$user_no_of_bases_to_increase_length)<length($species_seq))
									{
										$previous_spacer_bases_to_chop=substr($species_seq,$r_start-1-$user_no_of_bases_to_increase_length,$user_no_of_bases_to_increase_length); #--- the -1 is for substr(), its Zero based
									}
								else{
										#-------- make a dummy string and take substr of it : currently just ending the process----------------- Example: NC_019776
										
										if(($r_start-1) > 0)
											{
												$previous_spacer_bases_to_chop= substr($species_seq,0,$r_start-1);
												#print "**** $r_start-1 : $previous_spacer_bases_to_chop\n\n";
											}	
										
										
										#-------- make a dummy string and take substr of it : currently just ending the process----------------- Example: NC_019776
										#print "\n\n############## here #################\n\n";
										my $dummy_str="";
										for(my $i=1;$i<= ($user_no_of_bases_to_increase_length-length($previous_spacer_bases_to_chop));$i++)
											{
												$dummy_str=$dummy_str."-";
											}
										$previous_spacer_bases_to_chop=$dummy_str.$previous_spacer_bases_to_chop;
										#$previous_spacer_bases_to_chop=substr($species_seq,-$user_no_of_bases_to_increase_length);	##### its for CRISPRPriming project, where bases from the other end was taken
										
										
									}	
							}
						elsif($k1>4)
							{
								my $previous_line=$$modified_array[$#{$modified_array}]; chomp $previous_line; $previous_line=~s/\r+//g; $previous_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
						
								my @tmp_array2=split('\t',$previous_line);						
										
								$p_r_start=$tmp_array2[0];
								$p_r_seq=$tmp_array2[1]; my $p_r_seq_1=$p_r_seq;$p_r_seq_1=~s/-//g;
								
								if(defined $tmp_array2[2])
									{
										$p_s_seq=$tmp_array2[2];
									}
								else{
										$p_s_seq="-";
									}		
								
								my $p_s_seq_1=$p_s_seq;$p_s_seq_1=~s/-//g;
								
								if($tmp_array2[3])
									{
										$p_comment=$tmp_array2[3];
									}
								
								if($r_seq!~/^-/)
									{
										if(length($p_s_seq_1)>=$user_no_of_bases_to_increase_length)
											{
												#$current_spacer_bases_to_chop=substr($s_seq_1,0,$user_no_of_bases_to_increase_length);
												$previous_spacer_bases_to_chop=substr($p_s_seq_1,-$user_no_of_bases_to_increase_length);	
											}
										else{
												#--- add the available bases to $previous_spacer_bases_to_chop, and then add - to make it equal length on the left
												$previous_spacer_bases_to_chop=$p_s_seq_1;
												for(my $i=1;$i<=$user_no_of_bases_to_increase_length-length($p_s_seq_1);$i++) #Remember, starting from 1
													{
														$previous_spacer_bases_to_chop="-".$previous_spacer_bases_to_chop;
													}
											}
											
										#$previous_spacer_bases_to_chop=substr($p_s_seq,-$user_no_of_bases_to_increase_length);								
										
										
										
										
										$p_s_seq=~s/$previous_spacer_bases_to_chop$//;	
																
										my $modified_p_rec_line="$p_r_start\t$p_r_seq\t$p_s_seq\t$p_comment";
										$$modified_array[$#{$modified_array}]=$modified_p_rec_line;
									}
								else{
										$previous_spacer_bases_to_chop="";
										for(my $i=0;$i<$user_no_of_bases_to_increase_length;$i++)
											{
												$previous_spacer_bases_to_chop=$previous_spacer_bases_to_chop."-";
											}
									}
							}
							
						
						
						my $modified_r_start;
						if($r_seq!~/^-/)
							{
								$modified_r_start=$r_start-$user_no_of_bases_to_increase_length;
							}
						else{
								$modified_r_start=$r_start;
							}
						if($modified_r_start<1){$modified_r_start=1;}
						
						my $old_r_seq=&change_dots_to_bases($r_seq,$model_repeat);						
						
						my $modified_repeat="";						
						if(defined $previous_spacer_bases_to_chop and $previous_spacer_bases_to_chop=~/\S+/)
							{
								$modified_repeat=$previous_spacer_bases_to_chop.$old_r_seq;
							}
						else{
								$modified_repeat=$old_r_seq;
							}		
						
						my $modified_spacer=$s_seq;				
						
						
						my $modified_c_rec_line="";						
						if($k1<$#{$current_array}-1)
							{
								#$modified_rec_line="$r_start\t$modified_repeat\t$modified_spacer\t$comment";
								$modified_c_rec_line="$modified_r_start\t$modified_repeat\t$s_seq\t$comment";								
							}
						else{
								#$modified_rec_line="$r_start\t$modified_repeat\t|\t$comment";
								$modified_c_rec_line="$modified_r_start\t$modified_repeat\t|\t$comment";
							}
							
						#my $modified_c_rec_line="$modified_r_start\t$modified_repeat\t$s_seq\t$comment";
						
						#print "D:$modified_c_rec_line\n";
						push(@{$modified_array},$modified_c_rec_line);
						
						$new_r_seq=$modified_repeat;
					}
					
					
					
			#------------- store the repeats in a hash ------------------------------	
				if($hash_of_repeats{$new_r_seq})
					{
						$hash_of_repeats{$new_r_seq}=$hash_of_repeats{$new_r_seq}+1;
					}
				else{
						$hash_of_repeats{$new_r_seq}=1;
					}
						
				push(@arr_repeats,$new_r_seq);		
			}
		
		
		
		
				
		#----- now check for highest occurrence score of the modified repeats, and then find the model repeat ---
		my @arr_of_model_repeats;
		my $new_model_repeat="";
		my $heighest_repeat_occurrence_score=0;
		foreach my $repeat(sort{$hash_of_repeats{$b}<=>$hash_of_repeats{$a}} keys %hash_of_repeats)
			{
				#print "$user_side_to_increase_length:\t$repeat\n";
				if(not $arr_of_model_repeats[0])
					{
						#push(@arr_of_model_repeats,$repeat);
						$arr_of_model_repeats[0]=$repeat;
						$heighest_repeat_occurrence_score=$hash_of_repeats{$repeat};
					}
				elsif($hash_of_repeats{$repeat}>=$heighest_repeat_occurrence_score)
					{
						push(@arr_of_model_repeats,$repeat);
					}

			}
		
		
		#--------- check if there is just one repeat in the array
		if($#arr_of_model_repeats==0)
			{
				$new_model_repeat=$arr_of_model_repeats[0];
			}
		elsif($#arr_of_model_repeats>0)
			{
				#----- now select the best model repeat by scoring the array
				my %tmp_hash1;
				foreach my $repeat(@arr_of_model_repeats)
					{
						my $array_degeneracy_score=&get_array_degeneracy_score_for_undotted_repeats($repeat,\@arr_repeats);
						
						$tmp_hash1{$repeat}=$array_degeneracy_score
					}
				
				foreach my $repeat(sort{$tmp_hash1{$b}<=>$tmp_hash1{$a}}keys %tmp_hash1)
					{
						#print "$repeat\t$tmp_hash1{$repeat}\n";
						
						$new_model_repeat=$repeat;
						last;
					}							
			}
		else{
				#$new_model_repeat=$model_repeat;
						#$new_model_repeat=$arr_of_model_repeats[0];
			}
		
		
			

		#if($model_repeat ne $new_model_repeat)
		#	{
		#		$case_found=1;				
		
		
		
		#---- now re-create the array with the new model_repeat ------ 
				
		for(my $k1=0;$k1<=$#{$modified_array};$k1++)
			{
						#print "@{$current_array-[0]}\n";
						my $current_line=$$modified_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
									
						#print "A: $current_line\n";
									
						my @tmp_array=split('\t',$current_line);
						my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
								
						$r_start=$tmp_array[0];
						$r_seq=$tmp_array[1];
						$s_seq=$tmp_array[2];
						
						#---- special case: Fix added on 21-08-2013 
						if(not defined $s_seq or $s_seq eq "")
							{
								if($k1>0)
									{
										my $previous_line=$$modified_array[$k1-1]; chomp $previous_line; $previous_line=~s/\r+//g; $previous_line=~s/^\s+//;
										my @tmp_array2=split('\t',$previous_line);
										my $pre_s_seq=$tmp_array2[2];
										
										for(0..(length($pre_s_seq)-1))
											{
												$s_seq=$s_seq."-";
											}
										
									}
								else{	
										$s_seq="-";
									}	
							}
						#--------------------------------------------	
							
						if($tmp_array[3])
							{
								$comment=$tmp_array[3];
							}
						else{
								$comment="";
							}
						
						#my $old_r_seq=&change_dots_to_bases($r_seq,$model_repeat);
						my $new_r_seq=&change_bases_to_dots($r_seq,$new_model_repeat);
						
						$$modified_array[$k1]="$r_start\t$new_r_seq\t$s_seq\t$comment";	
						
						#print "D:$r_start\t$new_r_seq\t$s_seq\t$comment\n";	
			}
				
				#print "\nNew Consensus=$new_model_repeat\n\$model_repeat=$model_repeat\n";
				
				
				$model_repeat=$new_model_repeat;
		#	}		
		$case_found=1;
		return($model_repeat,$case_found);
	}



sub auto_detect_repeat_sides_and_bases_to_increase()
	{
		my($range,$repeat_extension_identity,$accession,$model_repeat,$current_array,$modified_array)=@_;
		
		my $case_found=0;
		my $user_side_to_increase_length="";
		my $alt_repeat_extension_identity=50;

				#print "Going to increase the repeat_lengths for an the array of $accession: with $user_side_to_increase_length,$user_no_of_bases_to_increase_length,$accession,$model_repeat\n\n";
		
		#---------- now open the sequence file and get the sequence string -------------
		open(SEQ,"$tmp_dir\/$accession\.fna") or print "$!";
		my @arr_seq=<SEQ>;
		close(SEQ);
		my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;
		#-------------------------------------------------------------------------------
		
		
		
		#------ first get max spacer length ------------------------------------------
		my $no_of_repeats=0;
		my $max_spacer_length=0;
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
				#print "@{$current_array-[0]}\n";
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
							
				#print "C: $current_line\n";
							
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq);
				my $comment="";
				my $current_spacer_bases_to_chop="";
						
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1]; my $r_seq_1=$r_seq;$r_seq_1=~s/-//g;
				
				if(not defined $tmp_array[2]){$tmp_array[2]="";}
				$s_seq=$tmp_array[2]; my $s_seq_1=$s_seq;$s_seq_1=~s/-//g;	#if(not defined $s_seq_1){print "\nError with $accession \n\n";}
				
				if(length($s_seq_1)>$max_spacer_length){$max_spacer_length=length($s_seq_1);}
				
				$no_of_repeats++;
			}
			
			
		#-------- now calculate the $repeat_extension_identity if it is set to AUTO
		if($repeat_extension_identity=~/AUTO/ or $repeat_extension_identity=~/^0/)
			{
				if($no_of_repeats<=3)
					{
						$repeat_extension_identity=100;	
						$alt_repeat_extension_identity=50;
					}
				elsif($no_of_repeats==4)
					{
						$repeat_extension_identity=75;
						$alt_repeat_extension_identity=50;
					}
				elsif($no_of_repeats==5)
					{
						$repeat_extension_identity=75;
						$alt_repeat_extension_identity=50;
					}
				elsif($no_of_repeats==6)
					{
						$repeat_extension_identity=80;
						$alt_repeat_extension_identity=50;
					}	
				else{
						$repeat_extension_identity=67;
						$alt_repeat_extension_identity=40;
					}			
			}	
			
		#------------------ get all the spacers in @arr_spacers
		my @arr_spacers;
		my @arr_spacers_rev;
		
		my %hash_of_repeats;
		my @arr_repeats;
		my $left_flanking_bases="";
		my $right_flanking_bases="";		
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
				#print "@{$current_array-[0]}\n";
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
							
				#print "C: $current_line\n";
							
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq);
				my $comment="";
				my $current_spacer_bases_to_chop="";
						
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1]; my $r_seq_1=$r_seq;$r_seq_1=~s/-//g;
				
				if(not defined $tmp_array[2]){$tmp_array[2]="";}
				$s_seq=$tmp_array[2]; my $s_seq_1=$s_seq;$s_seq_1=~s/-//g;	#if(not defined $s_seq_1){print "\nError with $accession \n\n";}
				
				
				
				if($k1==4)
					{
						
						push(@arr_spacers,$s_seq_1);
						
						my $t_s_seq_1=reverse($s_seq_1);
						push(@arr_spacers_rev,$t_s_seq_1);   			# or else current spacer will be missed
						
						if($no_of_repeats>2)#----- do not push the left flank if $no_of_repeats >2; i.e, at least 2 spacers are there to check
							{
								#--- get the left flank, and check if the end has bases or Ns, if Ns at the end then skip
								$left_flanking_bases=substr($species_seq,0,$r_start-1);
								
								if($left_flanking_bases=~/N$/i)
									{
										next;
									}
								
								#next;   
							}	
						
						
						
						if(($r_start-1-$max_spacer_length)>=0 and $r_start<length($species_seq))
							{
								
								#--- get the potential spacer from left flank: length equals to the first spacer or 100
								$left_flanking_bases=substr($species_seq,$r_start-1-$max_spacer_length,$max_spacer_length);
								#system("echo '$left_flanking_bases' >log.txt");
								
								if(not defined $left_flanking_bases){$left_flanking_bases="-";}
								
								$left_flanking_bases=reverse($left_flanking_bases);
								$s_seq_1=reverse($s_seq_1);
								
								unshift(@arr_spacers_rev,$left_flanking_bases);
							}
						else{
								
								#--- get the remaining sequences from left flank: length equals to the first spacer or 100
								$left_flanking_bases=substr($species_seq,0,$r_start-1);
								#system("echo '$left_flanking_bases' >log.txt");
								
								if(not defined $left_flanking_bases or $left_flanking_bases!~/\S/){$left_flanking_bases="-";}
								
								$left_flanking_bases=reverse($left_flanking_bases);
								$s_seq_1=reverse($s_seq_1);
								
								unshift(@arr_spacers_rev,$left_flanking_bases);
								
							}	
						
						
						
						
					}
				elsif($k1>4 and $k1<$#{$current_array}-1)
					{
						push(@arr_spacers,$s_seq_1);	
						
						#system("echo '$s_seq_1' >>log.txt");
						
						#-- now reverse the spacer bases and push in the array	
						
						$s_seq_1=reverse($s_seq_1);					
						push(@arr_spacers_rev,$s_seq_1);
					}
				elsif($k1==$#{$current_array}-1)
					{
						
						my $right_flanking_bases="";
						
						if(($r_start+length($r_seq_1)-1)<length($species_seq))
							{
								#--- get the right flank,
								if(length($species_seq)>=($r_start+length($r_seq_1)-1+$max_spacer_length) )
									{
										$right_flanking_bases=substr($species_seq,$r_start+length($r_seq_1)-1,$max_spacer_length);
									}
								else{	
										$right_flanking_bases=substr($species_seq,$r_start+length($r_seq_1)-1);
									}
								
								if($no_of_repeats>2)#----- do not push the right flank if $no_of_repeats >2; i.e, at least 2 spacers are there to check
									{
										#---  check if the begining has bases or Ns, if Ns at the start then skip							
											
										if($right_flanking_bases!~/\S+/ or $right_flanking_bases=~/^N/i)
											{
												next;
											}
										#push(@arr_spacers,$right_flanking_bases);	
										#next;   
									}
								else{	
									
										#--- get the potential spacer from right flank
										push(@arr_spacers,$right_flanking_bases);
																				
									}
							}
					}				
				
			}
		
		
		
		
		#my $alt_repeat_extension_identity=50+abs($#arr_spacers -2)*2; # so that with more number of rows, the minimum identity also increases		
		#print "\$repeat_extension_identity=$repeat_extension_identity\t\$alt_repeat_extension_identity=$alt_repeat_extension_identity\n";
		
		#----- get aligned spacers ----------------------------------------------------------
		#my @arr_aligned_spacers;	
		#&check_spacer_identity($accession,\@arr_spacers,\@arr_aligned_spacers);
		
		#my @arr_aligned_spacers_rev;	
		#&check_spacer_identity($accession,\@arr_spacers_rev,\@arr_aligned_spacers_rev);
		#exit;
		

		
		#----- get max spacer length and split the spacer bases in a 2D hash  ---------------		
		my %hash_of_spacer_bases;
		my $spacer_count=0;
		foreach my $spacer(@arr_spacers)
		#foreach my $spacer(@arr_aligned_spacers)
			{
				if(not defined $spacer){$spacer="";}
				
				chomp $spacer; $spacer=~s/\s+$//;
				#system("echo '$spacer_count $spacer' >>log.txt");
				#print "$spacer_count $spacer\n";
				my @arr_t1=split('',$spacer);
				my $base_count=0;
				foreach my $base(@arr_t1)
					{
						$hash_of_spacer_bases{$spacer_count}{$base_count}=$base;
						$base_count++;
					}
				
				if(length($spacer)>$max_spacer_length and $spacer_count<$#arr_spacers){$max_spacer_length=length($spacer);}		#------ dont include the final spacer (right flank) while getting the max spacer length
				
				$spacer_count++;
			}
		#print "\n";	
		
		#--- similarly do for the reversed ---
		my %hash_of_spacer_bases_rev;
		my $spacer_count1=0;
		foreach my $spacer(@arr_spacers_rev)
		#foreach my $spacer(@arr_aligned_spacers_rev)
			{
				chomp $spacer; $spacer=~s/\s+$//;
				#print "\tL: $spacer\n";
				#system("echo '$spacer_count1 $spacer' >>log.txt");
				my @arr_t1=split('',$spacer);
				my $base_count1=0;
				foreach my $base(@arr_t1)
					{
						$hash_of_spacer_bases_rev{$spacer_count1}{$base_count1}=$base;
						$base_count1++;
					}
				
				#if(length($spacer)>$max_spacer_length){$max_spacer_length=length($spacer);}		
				
				$spacer_count1++;
			}




					
		#------ now first check on the right side how far we can increase the repeat --------		
				#system("echo '$repeat_extension_identity' >log.txt");		
		my $right_extension=0;		
		my $right_stop_found=0;
		for(my $i=0;$i<$max_spacer_length;$i++)
			{
				#print "$i:\n";	
				#--- make a string of all the bases in the current and next column
				my $current_column_str="";
				my $next_column_str="";
				my $next_column_str2="";		
				for(my $j=0;$j<$spacer_count;$j++)
					{
						my $base="-";
						my $next_base="-";
						my $next_base2="-";
						
						#---- current column
						if(defined $hash_of_spacer_bases{$j}{$i})
							{
								$base=$hash_of_spacer_bases{$j}{$i};
							}						
						$current_column_str=$current_column_str.$base;
						
						#----- next column
						if(defined $hash_of_spacer_bases{$j}{$i+1})
							{
								$next_base=$hash_of_spacer_bases{$j}{$i+1};
							}						 
						$next_column_str=$next_column_str.$next_base;
						#----- second next column --------------------------------		
						if(defined $hash_of_spacer_bases{$j}{$i+2})
							{
								$next_base2=$hash_of_spacer_bases{$j}{$i+2};
							}						 
						$next_column_str2=$next_column_str2.$next_base2;	
						#print "\t$j:\$next_column_str2=$next_column_str2\n";
								
					}
					
	
				my($the_base_with_highest_occurrence_cur,$occurrence_score_current_base)=&get_the_base_with_highest_occurrence_from_a_string($current_column_str,$spacer_count);				
				my($the_base_with_highest_occurrence,$occurrence_score_next_base)=&get_the_base_with_highest_occurrence_from_a_string($next_column_str,$spacer_count);
				my($the_base_with_highest_occurrence2,$occurrence_score_next_base2)=&get_the_base_with_highest_occurrence_from_a_string($next_column_str2,$spacer_count);
				#--- now check if the base with highest occurrence is >= $repeat_extension_identity
				
				if($the_base_with_highest_occurrence_cur eq "-"){$right_stop_found=1;last;}
				
				if( ($occurrence_score_current_base >= $repeat_extension_identity) or ($spacer_count >=3  and $spacer_count1 <5 and length($right_flanking_bases)==1 and $occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>=$alt_repeat_extension_identity and $occurrence_score_next_base2>=$alt_repeat_extension_identity) or ($spacer_count >=3 and $occurrence_score_current_base >$alt_repeat_extension_identity and $right_flanking_bases =~/-/) or ($occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>=$repeat_extension_identity and $occurrence_score_next_base2>=$repeat_extension_identity) or ($spacer_count >=3 and $occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>$alt_repeat_extension_identity and $occurrence_score_next_base2>$alt_repeat_extension_identity))
					{														
						$right_extension++;
					}
				else{
						#print "Right Base index:$i \$occurrence_score_current_base=$occurrence_score_current_base \$repeat_extension_identity=$repeat_extension_identity \t $the_base_with_highest_occurrence_cur,$occurrence_score_current_base $the_base_with_highest_occurrence,$occurrence_score_next_base $the_base_with_highest_occurrence2,$occurrence_score_next_base2\n";
						$right_stop_found=1;
					}
				
						
				if($right_stop_found==1){last;}		
			}
		
		if($right_extension>0){$user_side_to_increase_length="RIGHT-$right_extension,";}
		#print "\$user_side_to_increase_length=$user_side_to_increase_length\n";
		
		
		
		
		#------ now check on the left side to see how far we can increase the repeat -------- 
		
		#system("echo '$repeat_extension_identity' >log.txt");
		
		my $left_extension=0;
		my $left_stop_found=0;
		for(my $i=0;$i<$max_spacer_length-$right_extension;$i++)
			{

				
				#--- make a string of all the bases in the current and next column
				my $current_column_str="";
				my $next_column_str="";
				my $next_column_str2="";
				for(my $j=0;$j<$spacer_count1;$j++)
					{
						my $base="-";
						my $next_base="-";
						my $next_base2="-";
						
						if(defined $hash_of_spacer_bases_rev{$j}{$i})
							{
								$base=$hash_of_spacer_bases_rev{$j}{$i};
							}						
						$current_column_str=$current_column_str.$base;
						
						#----- next column ---------------------------------------
						if(defined $hash_of_spacer_bases_rev{$j}{$i+1})
							{
								$next_base=$hash_of_spacer_bases_rev{$j}{$i+1};								
							}						 
						$next_column_str=$next_column_str.$next_base;
						#----- second next column --------------------------------		
						if(defined $hash_of_spacer_bases_rev{$j}{$i+2})
							{
								$next_base2=$hash_of_spacer_bases_rev{$j}{$i+2}
							}						
						$next_column_str2=$next_column_str2.$next_base2;		
					}
				#--- now check if the base with highest occurrence is >= $repeat_extension_identity
				my($the_base_with_highest_occurrence_cur,$occurrence_score_current_base)=&get_the_base_with_highest_occurrence_from_a_string($current_column_str,$spacer_count1);				
				my($the_base_with_highest_occurrence,$occurrence_score_next_base)=&get_the_base_with_highest_occurrence_from_a_string($next_column_str,$spacer_count1);
				my($the_base_with_highest_occurrence2,$occurrence_score_next_base2)=&get_the_base_with_highest_occurrence_from_a_string($next_column_str2,$spacer_count1);
				#print "\n******\n($the_base_with_highest_occurrence2,$occurrence_score_next_base2)=get_the_base_with_highest_occurrence_from_a_string($next_column_str2,$spacer_count1)\n";
				#--- now check if the base with highest occurrence is >= $repeat_extension_identity
				#print qq~ ($occurrence_score_current_base >= $repeat_extension_identity) or \
				#	($spacer_count1 >=3 and length($left_flanking_bases)==1 and $occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>=$alt_repeat_extension_identity and $occurrence_score_next_base2>=$alt_repeat_extension_identity) or \
				#	($spacer_count1 ==3 and $occurrence_score_current_base >$alt_repeat_extension_identity and $left_flanking_bases =\~/-/) or \
				#	($occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>=$repeat_extension_identity and $occurrence_score_next_base2>=$repeat_extension_identity)
				#	~;
				
				if($the_base_with_highest_occurrence_cur eq "-"){$left_stop_found=1;last;}
				
				#if($spacer_count1 >=3 and $spacer_count1 <5  and length($left_flanking_bases)==1 and $occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>=$alt_repeat_extension_identity and $occurrence_score_next_base2>=$alt_repeat_extension_identity)
				#	{
				#		print "A: This is the error: $i\n $spacer_count1 >=3 and length($left_flanking_bases)==1 and $occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>=$alt_repeat_extension_identity and $occurrence_score_next_base2>=$alt_repeat_extension_identity\n";
				#	}
				#if($occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>=$repeat_extension_identity and $occurrence_score_next_base2>=$repeat_extension_identity)
				#	{
				#		print "B: This is the error: $i\n $occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>=$repeat_extension_identity and $occurrence_score_next_base2>=$repeat_extension_identity\n";
				#	}
				
				#print "if(($occurrence_score_current_base >= $repeat_extension_identity) or ($spacer_count1 >=3 and $spacer_count1 <5 and length($left_flanking_bases)==1 and $occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>=$alt_repeat_extension_identity and $occurrence_score_next_base2>=$alt_repeat_extension_identity) or ($spacer_count1 ==3 and $occurrence_score_current_base >$alt_repeat_extension_identity and $left_flanking_bases =~/-/) or ($occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>=$repeat_extension_identity and $occurrence_score_next_base2>=$repeat_extension_identity))\n";	
				#if(($occurrence_score_current_base >= $repeat_extension_identity) or ($spacer_count1 >=3 and length($left_flanking_bases)==1 and $occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>=$alt_repeat_extension_identity and $occurrence_score_next_base2>=$alt_repeat_extension_identity) or ($spacer_count1 ==3 and $occurrence_score_current_base >$alt_repeat_extension_identity and $left_flanking_bases =~/-/) or ($occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>=$repeat_extension_identity and $occurrence_score_next_base2>=$repeat_extension_identity))
				if(($occurrence_score_current_base >= $repeat_extension_identity) or ($spacer_count1 >=3 and $spacer_count1 <5 and length($left_flanking_bases)==1 and $occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>=$alt_repeat_extension_identity and $occurrence_score_next_base2>=$alt_repeat_extension_identity) or ($spacer_count1 ==3 and $occurrence_score_current_base >$alt_repeat_extension_identity and $left_flanking_bases =~/-/) or ($occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>=$repeat_extension_identity and $occurrence_score_next_base2>=$repeat_extension_identity) or ($spacer_count >=3 and $occurrence_score_current_base >=$alt_repeat_extension_identity and $occurrence_score_next_base>$alt_repeat_extension_identity and $occurrence_score_next_base2>$alt_repeat_extension_identity))
					{								
						$left_extension++;
					}
				else{
						$left_stop_found=1;
					}						
					
				if($left_stop_found==1){last;}		
			}
		
		if($left_extension>0){$user_side_to_increase_length="LEFT-$left_extension,".$user_side_to_increase_length;}
		
		if($user_side_to_increase_length eq "")
			{
				$user_side_to_increase_length="NA-0,";
			}	
		#system("echo '$user_side_to_increase_length $spacer_count $max_spacer_length' >>log.txt");		
		
		$case_found=1;
		return($user_side_to_increase_length,$case_found);
	}



sub check_spacer_identity()
	{
		my($accession,$arr_spacers,$arr_aligned_spacers)=@_;
		my $spacers_identity=0;
		
		#---- create a fasta file with only the spacers with some sequence----
		my $time_1 = &get_unique_id();
		
		my $tmp_file=		$accession.$time_1."_tmp_spacers.txt";
		my $tmp_output=		$accession.$time_1."_tmp_spacers.aln";
		my $tmp_output_dnd=	$accession.$time_1."_tmp_spacers.dnd";
		
		open(FA,">$tmp_dir\/$tmp_file") or print "$!";
		flock(FA,2);
		
		for(my $i=0;$i<=$#{$arr_spacers};$i++)
			{
				my $t_spacer=$$arr_spacers[$i];
				if($t_spacer eq "")
					{
						$t_spacer="-";
					}
						
				print FA ">S_$i\n$t_spacer\n";				
			}
		close(FA);
		
		#---- run clustalW to get the repeat and spacers
		
		
		system("clustalw -INFILE=$tmp_dir\/$tmp_file -OUTFILE=$tmp_dir\/$tmp_output -QUIET -ALIGN -NUMITER=50 -ENDGAPS -OUTORDER=INPUT >/dev/null 2>&1");
		
		
		open(RD,"$tmp_dir\/$tmp_output");
		my @arr_output=<RD>;
		close(RD);
		
		
		unlink("$tmp_dir\/$tmp_file");	
		unlink("$tmp_dir\/$tmp_output");
		unlink("$tmp_dir\/$tmp_output_dnd");
		#--------------------------------------------------------------------------
		my %hash_of_lines;
		#foreach my $line(@arr_output)
		
		for(my $i=0;$i<=$#arr_output;$i++)
			{
				my $line=$arr_output[$i];
				chomp $line; $line=~s/\r$//;
				
				if($line=~/^S_/)
					{	
											
						
						$line=~s/^S_//; $line=~s/\s+/\t/g;
						my @arr_t1=split('\t',$line);
						
						if(defined $hash_of_lines{$arr_t1[0]})
							{
								$hash_of_lines{$arr_t1[0]}=$hash_of_lines{$arr_t1[0]}.$arr_t1[1];
							}
						else{
								
								$hash_of_lines{$arr_t1[0]}=$arr_t1[1];
							}	
					}
				#elsif($line=~/^\s{16}/)
				#	{
				#		$line=~s/^\s{16}//;
				#		
				#		my $hash_index=$#arr_spacers+1;
				#		if($arr_output[$i-1]=~/^S_(\d+)/){$hash_index=$1+1;}
				#		
				#		if(defined $hash_of_lines{$hash_index})
				#			{
				#				$hash_of_lines{$hash_index}=$hash_of_lines{$hash_index}.$line;
				#			}
				#		else{								
				#				if($line ne "")
				#					{
				#						$hash_of_lines{$hash_index}=$line;
				#					}	
				#			}
				#	}	
					
			}
		
		my $starting_gap=0;
		my $model_repeat="";
		foreach my $seq_index(sort{$a<=>$b} keys %hash_of_lines)
			{
				#print "$seq_index\t$hash_of_lines{$seq_index}\n";
				
				$$arr_aligned_spacers[$seq_index]=$hash_of_lines{$seq_index};
				
			}
			

		
		return 1;
	}



sub get_the_base_with_highest_occurrence_from_a_string()
	{
		my ($input_string,$spacer_count)=@_;
		my @t_arr=split('',$input_string);
		my %t_hash;
		foreach my $base(@t_arr)
			{
				if(defined $t_hash{$base})
					{
						$t_hash{$base}=$t_hash{$base}+1;
					}
				else{
						$t_hash{$base}=1;
					}
			}
							
		my $the_base_with_highest_occurrence="-";
		my $occurrence_score=0;
							
		foreach my $base(sort{$t_hash{$b}<=>$t_hash{$a}} keys %t_hash)
			{
				my $occurrence=$t_hash{$base};						
				my $occurrence_p=sprintf("%.0f",($occurrence*100)/$spacer_count);
								
				$the_base_with_highest_occurrence=$base;
				$occurrence_score=$occurrence_p;
				last;
			}
		return($the_base_with_highest_occurrence,$occurrence_score);		
	}


sub obtain_array_scores()
	{
		my($current_repeat,$model_repeat,$div_length)=@_;
		
		my @arr_1=split('',$current_repeat);
		my @arr_2=split('',$model_repeat);
		
		#my $array_score=0;
		my $left_degeneracy_score=0;
		my $center_degeneracy_score=0;
		my $right_degeneracy_score=0;
		
		for(my $i=0;$i<=$#arr_2;$i++)   
			{
				if($arr_1[$i] ne ".")
					{						
						if($i>=0 and $i<=$div_length)
							{
								if($i==0){$left_degeneracy_score--;next;}
								
								
								if($arr_1[$i] eq "-" and $arr_1[$i-1] eq "-")
									{
										#$left_degeneracy_score--;
									}
								else{									
										$left_degeneracy_score--;
									}	
									
							}
						elsif($i>$div_length and $i<=(2*$div_length))
							{
								#if($arr_1[$i] ne "." and $arr_1[$i-1] ne "."){next;}
								#else{
										$center_degeneracy_score--;
								#	}
							}
						else{								
								if($arr_1[$i] eq "-" and $arr_1[$i-1] eq "-")
									{
										#$right_degeneracy_score--;
									}
								else{
										$right_degeneracy_score--;
									}	
									
							}
						#-----	
						#$array_score=$array_score--;			
					}
			}
		
		return ($left_degeneracy_score,$center_degeneracy_score,$right_degeneracy_score);
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


sub change_dots_to_bases()
	{
		my($r_string,$mr_string)=@_;
		
		if($r_string eq "" or $mr_string eq "" ){return("");}
		#print "\$bf_string,\$mr_string $bf_string,$mr_string\n";
		my $return_string;
		
		my @arr_1=split('',$r_string);
		my @arr_2=split('',$mr_string);
		
		for(my $i=0;$i<=$#arr_2;$i++)
			{
				if($arr_1[$i] eq ".")
					{
						$return_string=$return_string.$arr_2[$i];
					}
				else{
						$return_string=$return_string.$arr_1[$i];
					}	
			}
		
		return($return_string);
	}



sub search_unidentified_repeat_in_spacers()
	{		
		
		my($range,$su_dynamic_search,$allowed_percent_similarity,$accession,$model_repeat,$avg_spacer_length,$current_array,$modified_array)=@_;
		
		#print "Going to search for unidentified repeat in spacer sequences for $accession:\n\n";
		
		my $case_found=0;
		my $left_flank="";
		my $right_flank="";

		

		my $devision_length=int(length($model_repeat)/4);
		
		my $model_spacer_length=28;
		
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
				
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;							
						#print "A: $current_line\n";
							
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
						
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1];
				$s_seq=$tmp_array[2];
				
				my $existing_insertion_bases="";
				my $existing_insertion_positions="";
				if($tmp_array[3])
					{
						$comment=$tmp_array[3]; $comment=~s/^\s+//;
						if($tmp_array[3]!~/^\s{0,1}Del/)
							{
								($existing_insertion_bases,$existing_insertion_positions)=split(' ',$comment);
								if($existing_insertion_positions)
									{
										$existing_insertion_positions=~s/\[//g;
										$existing_insertion_positions=~s/\]//g;
									}
							}	
					}
				else{
						$comment="";
					}	
						
																				
				$r_length=length($r_seq);
				
				#------------ check if there are no base(s) -----
				if(not defined $s_seq or $s_seq eq ""){$s_seq="-";}							
				$s_length=length($s_seq);
				

				my $gapless_s_seq=$s_seq; $gapless_s_seq=~s/-//g;
				if(not defined $gapless_s_seq or $gapless_s_seq eq "")
					{
						push(@{$modified_array},$current_line);
						next;
					}
				#------------------------------------------------
				#system("echo '$avg_spacer_length' >>log.txt");
				
				#print "\n\n\n\n$s_length\n\n\n";
				if((length($gapless_s_seq)>=$avg_spacer_length*1.2 and length($gapless_s_seq)>=length($model_repeat)) and ($k1<$#{$current_array}-1) and $k1>=4) #if($s_length>=($avg_spacer_length*1.5) and ($k1<$#{$current_array}-1) and $k1>=4)
					{
						push(@{$modified_array},$current_line);	
						
						#next;
						
						#print "\nInside: with $current_line\n\n\n";
						#$s_seq=$right_flank;
						
						my $ref_repeat_seq;
						my $ref_spacer_seq;
						
						
						if($k1>4)
							{
								my @arr_previous_rec=split('\t',$$modified_array[$#{$modified_array}-1]);		#split('\t',$$current_array[$k1-1]);												
								
								#my $last_repeat_start=$arr_previous_rec[0];
								$ref_repeat_seq=$arr_previous_rec[1];
								if(not defined $arr_previous_rec[2])
									{
										$ref_spacer_seq="";
									}
								else{	
										$ref_spacer_seq=$arr_previous_rec[2];  #--- remember, ref_spacer_seq will always be 0
									}	
							}
						else{
							
								
								my @arr_next_rec=split('\t',$$current_array[$k1+1]);		#split('\t',$$current_array[$k1-1]);												
								
										#my $last_repeat_start=$arr_previous_rec[0];
								$ref_repeat_seq=$arr_next_rec[1];
								
								if(not defined $arr_next_rec[2])
									{
										$ref_spacer_seq="";
									}
								else{	
										$ref_spacer_seq=$arr_next_rec[2];  #--- remember, ref_spacer_seq will always be 0
									}	
							}	
						
						#--- use while condition -----------------------------------------------------------
						my ($original_best_fitting_pattern,$best_fitting_pattern,$insertion_base,$position);
						
						$best_fitting_pattern="-";
						
						my $pos;
						#my $next_spacer_seq="";
						my $repeat_seq_before_flank=$r_seq;
						
						my $skip_this_rec=0;
						my $extension_found=0;
						my $potential_next_repeat="";
						my $potential_previous_spacer="";
						
						my $no_of_tries=0;
						while($best_fitting_pattern ne "")
							{
							
								$no_of_tries++;
								if($no_of_tries>5){last;}
								
								#-------- find the repeat_start position from the last record -----------
								my @arr_cur_rec=split('\t',$$modified_array[$#{$modified_array}]);	
								
								#print "\t$k1: @arr_cur_rec\n";
								my $cur_r_start=$arr_cur_rec[0];
								my $cur_r_seq=$arr_cur_rec[1];											
								my $cur_s_seq="";
								if($arr_cur_rec[2])
									{
										$cur_s_seq=$arr_cur_rec[2]; $cur_s_seq=~s/-//g;
									}	
								
								if($cur_s_seq eq ""){last;}
								
								my $cur_repeat_seq=&change_dots_to_bases($cur_r_seq,$model_repeat);
								
								#print "$model_repeat\n$cur_repeat_seq\n\n";
								
								#if($cur_s_seq !=0)
								#	{
								#		$skip_this_rec=1;
								#		last;
								#	}
								my $cur_comment=$arr_cur_rec[3];
								
								my $cur_no_of_insertions=0;
								my ($cur_insertion_bases,$cur_insertion_positions);
								if($arr_cur_rec[3])
									{
										#$comment=$tmp_array[3];
										my $c_comment=$cur_comment; $c_comment=~s/^\s+//;
										#if($c_comment=~/^Insertions/)
										#	{
										#		$c_comment=~ s/Insertions of //;
										#	}
										#els
										if($c_comment!~/^Del/)
											{
												#next;											
											
												my @tmp_array=split(' ',$c_comment);
												
												#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
												$cur_insertion_bases=$tmp_array[0];
												$cur_insertion_positions=$tmp_array[1];
												#print "\$cur_insertion_bases=$cur_insertion_bases\n";
												
												$cur_insertion_bases=~s/,//g;
												$cur_no_of_insertions=length($cur_insertion_bases);
												#print "\$cur_no_of_insertions=$cur_no_of_insertions\n";										
										}
									}							

								#--------------------------------------------------------------------------------------------------				
								if($ref_spacer_seq=~/^\|/ or $ref_spacer_seq=~/^\0/ or $ref_spacer_seq=~/^NNNNNNNNNNN/)
									{
										$ref_spacer_seq=$cur_s_seq;
									}
								my $string_length=int((length($ref_repeat_seq) + $model_spacer_length)*1.25);
								
								#if($string_length<=length($cur_s_seq))
								#	{								
								#		$potential_previous_spacer=substr($cur_s_seq,0,$string_length);
								#	}
								#else{
										$potential_previous_spacer=$cur_s_seq;
								#	}		
								#print "\$potential_previous_spacer=$potential_previous_spacer\n";
								
								
								
								my($l_spacer,$l_spacer_pos);
								
								my $reference_repeat=$model_repeat;
								if($su_dynamic_search==1 and ($cur_r_seq!~/-/ or not $cur_no_of_insertions >0))
									{
										$reference_repeat=&change_dots_to_bases($cur_r_seq,$model_repeat);
									}

								#----  static method: searched with model repeat -------------------------------
								($original_best_fitting_pattern,$best_fitting_pattern,$l_spacer,$l_spacer_pos,$insertion_base,$position)=&get_best_fitting_string_in_flanks($range,$accession,$allowed_percent_similarity,$potential_previous_spacer,$reference_repeat,"RIGHT",0);
								#print "\nR. $original_best_fitting_pattern,$best_fitting_pattern,$l_spacer,$l_spacer_pos,$insertion_base,$position\n";	
								#system("echo '$original_best_fitting_pattern,$best_fitting_pattern,$l_spacer,$l_spacer_pos,$insertion_base,$position' >>log.txt");
								
								if($l_spacer eq "-"){$l_spacer="";}
								
								
								my $similarity_score;
								
								if(not $best_fitting_pattern)
									{
										$best_fitting_pattern="";
										$similarity_score=0;
									}
								else{
										$similarity_score=&get_similarity_score($best_fitting_pattern,$model_repeat);
										#print "\$similarity_score=$similarity_score\n";
										#---- find out how many gaps are there ---
										my $no_gaps=$best_fitting_pattern=~s/-/-/g;
										if($no_gaps>1)
											{
												$similarity_score=$similarity_score -$no_gaps*2;
											}
										if($insertion_base)
											{
												my $no_of_ins= $insertion_base=~s/,/,/g;
												#print "\$insertion_base==$insertion_base\t $similarity_score=$similarity_score-$no_of_ins*2\t$similarity_score<(length($model_repeat)*$allowed_percent_similarity/100)\n";
												$similarity_score=$similarity_score-$no_of_ins*2;
											}
									}	
								
								
								#--------------------------------------------------------------------------------	
								#----- check best fitting pattern to have more than 60% similarity with $model repeat
								if($best_fitting_pattern eq "" or $similarity_score<(length($model_repeat)*$allowed_percent_similarity/100) or length($l_spacer)>($model_spacer_length*1.5))
									{	
										#print "Skipping...\n";									
										$skip_this_rec=1;#------ no need to push anything as the record is already there -----------------------
										last;
									}												
								
								
								
								#------------------------ now create the two lines which will be pushed
									
								#------------ increase the $extension_found -----------------------------------------
								$case_found=1;
								$extension_found++;	
								$potential_next_repeat=&change_bases_to_dots($best_fitting_pattern,$model_repeat);
								#------------------------------------------------------------------------------------
								#my $pos;
								my $current_repeat=$cur_r_seq;		
								my $next_repeat=&change_bases_to_dots($best_fitting_pattern,$model_repeat);
								
								
								#-------- now store the new rec ------------------------------------------------
								#my $gapless_p_r_seq=$potential_next_repeat; $gapless_p_r_seq=~s/-//g;
								my $gapless_p_r_seq=$cur_r_seq; $gapless_p_r_seq=~s/-//g;
								my $gapless_p_s_seq=$l_spacer; $gapless_p_s_seq=~s/-//g;						
								
								
								my $comma_removed_insertion_bases=$insertion_base; 	$comma_removed_insertion_bases=~s/,//g;
									
								my $next_repeat_start=$cur_r_start+length($gapless_p_r_seq)+length($gapless_p_s_seq)+$cur_no_of_insertions;#+length($comma_removed_insertion_bases);	
								
								#print "\$next_repeat_start [$next_repeat_start] = $cur_r_start + length($gapless_p_r_seq) + length($gapless_p_s_seq) + $cur_no_of_insertions\n";
								#my $current_spacer=0;
								
								#print "OK\n";
								
								#--------------------------------------------------------------------------								
								my $comment1="";						
								if(length($l_spacer)<($model_spacer_length*0.80))
									{	
										if($cur_insertion_bases)
											{							
												$comment1=$cur_comment." Deletion [$next_repeat_start]";													
											}
										else{
												$comment1="Deletion [$next_repeat_start]";	
											}		
									}
								else{
										if($cur_insertion_bases)
											{							
												$comment1=$cur_comment;													
											}
										else{
												$comment1="";	
											}
									}	
								
								#------------ if l_spacer is blank, fill it with dashes -----------------
								if(not $l_spacer)
									{
										$l_spacer="";										
										$pos=$r_start+$r_length;
										
										#print "OK :",length($cur_s_seq),"\n";
										
										for(my $i=0;$i<$model_spacer_length;$i++)
											{											
												$l_spacer=$l_spacer."-";
											}
									}
								else{
										$pos=$next_repeat_start;		#-- deletion always happens just befor the next_start								
										$cur_s_seq=$l_spacer;
									}
									
								$cur_s_seq=~s/-//g;		
													
								#------- modify the previous record -------------------------------------
									
								my $modified_previous_rec="$cur_r_start\t$cur_r_seq\t$l_spacer\t$comment1";	
								$$modified_array[$#{$modified_array}]=$modified_previous_rec;
								
								#print "$modified_previous_rec\n";
								
								
								#-------------------------------------------------------------------------  
								#-------- now get the actual base postions for the insertions -------------
								my @array_positions=split(',',$position);
								my $new_position="";
								foreach my $p(@array_positions)
									{
										my $p1=$next_repeat_start+$p;
										$new_position=$new_position.",".$p1;
									}
								$new_position=~s/^,//;
								$position=$new_position;  
								my $comment2=""; 							
								if($insertion_base ne "")
									{
										$comment2="$insertion_base [$position]";	
									}
									
									
									
								#---------------------------------- now shorten the s_seq----------------------------------------
								my $gapless_best_fitting_pattern=$original_best_fitting_pattern;
								   $gapless_best_fitting_pattern=~s/-//g;
								   
								my $chopping_pattern= $l_spacer.$gapless_best_fitting_pattern; 
								$chopping_pattern=~s/-//g;
								
								#print "\$chopping_pattern$chopping_pattern\n";
								
								$s_seq=~s/^$chopping_pattern//;	
									
								#if($s_seq eq ""){$s_seq="-";}
								if($s_seq eq "")
									{
										$s_seq="";									
										
										#print "OK :",length($cur_s_seq),"\n";
										
										for(my $i=0;$i<$model_spacer_length;$i++)
											{											
												$s_seq=$s_seq."-";
											}
									}
								$s_seq=~s/-+$/-/;						
								my $new_rec_line1="$next_repeat_start\t$potential_next_repeat\t$s_seq\t$comment2";
								#print "\$new_rec_line1=$new_rec_line1\n\n";
								
								push(@{$modified_array},$new_rec_line1);
								
								
								
								#$s_seq=$right_flank;
								
								#exit;							
							} #--- end of while
						
						if($skip_this_rec==1){next;}

					}
				else{
						#my $new_rec_line="$r_start\t$r_seq\t$s_seq\t$comment";
						push(@{$modified_array},$current_line);
					}										
			}
		return($case_found);		
	}




sub get_similarity_score()
	{
		my($current_repeat,$model_repeat)=@_;
		
		
		
		my $similarity_score=0;
		
		if(not defined $current_repeat or $current_repeat eq "" or $model_repeat eq ""){return $similarity_score;}
		
		
		my @arr_1=split('',$current_repeat);
		my @arr_2=split('',$model_repeat);
		
		
		
		for(my $i=0;$i<=$#arr_2;$i++)   
			{
				if(defined $arr_1[$i] and $arr_1[$i] eq $arr_2[$i])
					{
						$similarity_score=$similarity_score+1;	
					}	
			}
		
		return $similarity_score;
	}






sub check_consensus_sequence()
	{
		my($range,$accession,$model_repeat,$avg_spacer_length,$current_array,$modified_array)=@_;
		#print "Going to check for consensus sequence with $model_repeat=$model_repeat in $accession:\n\n";
		
		my $case_found=0;
		my $new_model_repeat="";
		
		
		#---------- now open the sequence file and get the sequence string -------------
		#open(SEQ,"$tmp_dir\/$accession\.fna") or print "$!";
		#my @arr_seq=<SEQ>;
		#close(SEQ);
		#my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;
		
		my @arr_repeats;
		my %hash_of_repeats;
		
		my @arr_model_repeat_bases=split('',$model_repeat);
		
		#print "\n";
				
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
				#print "@{$current_array-[0]}\n";
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
							
				#print "C: $current_line\n";
							
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
						
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1];
				
				if($tmp_array[2])
					{
						$s_seq=$tmp_array[2];
					}
				else{
						$s_seq="";
					}		
				
				if($tmp_array[3])
					{
						$comment=$tmp_array[3];
					}
				else{
						$comment="";
					}
																				
				$r_length=length($r_seq);
				$s_length=length($s_seq);
				
				
				my $new_r_seq=&change_dots_to_bases($r_seq,$model_repeat);
				
				if(defined $hash_of_repeats{$new_r_seq})
					{
						$hash_of_repeats{$new_r_seq}=$hash_of_repeats{$new_r_seq}+1;
					}
				else{
						$hash_of_repeats{$new_r_seq}=1;
					}
					
					
				push(@arr_repeats,$new_r_seq);	
				
					
			}
		
		my @arr_of_model_repeats;
		my $heighest_repeat_occurrence_score=0;
		foreach my $repeat(sort{$hash_of_repeats{$b}<=>$hash_of_repeats{$a}} keys %hash_of_repeats)
			{
				if(not $arr_of_model_repeats[0])
					{
						#push(@arr_of_model_repeats,$repeat);
						$arr_of_model_repeats[0]=$repeat;
						$heighest_repeat_occurrence_score=$hash_of_repeats{$repeat};
					}
				elsif($hash_of_repeats{$repeat}>=$heighest_repeat_occurrence_score-0)
					{
						push(@arr_of_model_repeats,$repeat);
						#unshift(@arr_of_model_repeats,$repeat);
						
						#$arr_of_model_repeats[0]=$repeat;
					}
	
				#print "$repeat\t$hash_of_repeats{$repeat}\n";		
				#$new_model_repeat=$repeat;
				
				#print "\t\t\$new_model_repeat=$new_model_repeat\n";
				#last;
			}
		
		
		#--------- check if there is just one repeat in the array
		if($#arr_of_model_repeats==0)
			{
				$new_model_repeat=$arr_of_model_repeats[0];
				#print "Only one repeat got pushed\n";
			}
		elsif($#arr_of_model_repeats>0)
			{
				#----- now select the best model repeat by scoring the array
				my %tmp_hash1;
				foreach my $repeat(@arr_of_model_repeats)
					{
						my $array_degeneracy_score=&get_array_degeneracy_score($model_repeat,$repeat,\@arr_repeats);
						
						$tmp_hash1{$repeat}=$array_degeneracy_score;
						#system("echo '$repeat $array_degeneracy_score' >>log1.txt");
						
						#print "\$tmp_hash1{$repeat}=$array_degeneracy_score\n\n";
					}
				
				foreach my $repeat(sort{$tmp_hash1{$b}<=>$tmp_hash1{$a}}keys %tmp_hash1)
					{
						#print "$repeat\t$tmp_hash1{$repeat}\n";
						if($repeat!~/-/)
							{
								$new_model_repeat=$repeat;
								last;
							}	
					}	
					
									
			}
		else{
				$new_model_repeat=$model_repeat;
			}
		
		
			

		
		
		#if($new_model_repeat eq "" or uc($model_repeat) eq uc($new_model_repeat) or length($model_repeat)!=length($new_model_repeat)){return($left_flank,$right_flank,$model_repeat,$case_found);}
		
		if($model_repeat ne $new_model_repeat and $new_model_repeat ne "")
			{
				$case_found=1;				
				#---- now re-create the array with the new model_repeat ------ 
				
				for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
					{
						#print "@{$current_array-[0]}\n";
						my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
									
						#print "A: $current_line\n";
									
						my @tmp_array=split('\t',$current_line);
						my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
								
						$r_start=$tmp_array[0];
						$r_seq=$tmp_array[1];
						
						if(defined $tmp_array[2])
							{
								$s_seq=$tmp_array[2];
							}
						else{
								$s_seq="-";
							}
							
						if($tmp_array[3])
							{
								$comment=$tmp_array[3];
							}
						else{
								$comment="";
							}
						
						my $old_r_seq=&change_dots_to_bases($r_seq,$model_repeat);
						my $new_r_seq=&change_bases_to_dots($old_r_seq,$new_model_repeat);
						
						$$current_array[$k1]="$r_start\t$new_r_seq\t$s_seq\t$comment";	
						
						#print "D:$r_start\t$new_r_seq\t$s_seq\t$comment\n";	
					}
				
				#print "New Consensus=$new_model_repeat\n\$model_repeat=$model_repeat\n";
				
				
				$model_repeat=$new_model_repeat;
			}
			
			
			#print "\$model_repeat=$model_repeat\n\n";
				
			return($model_repeat,$case_found);		
	}


sub get_array_degeneracy_score()
	{
		my($model_repeat,$current_repeat,$ref_arr_repeats)=@_;
		
		my $arr_degen_score=0;
		
		foreach my $repeat(@{$ref_arr_repeats})
			{
				my $new_r_seq=&change_dots_to_bases($repeat,$model_repeat);
				
				my $sim_score=&get_similarity_score($new_r_seq,$current_repeat);
				
				my $neg_score=length($current_repeat)-$sim_score;
				
				$arr_degen_score=$arr_degen_score-$neg_score;
			}
			
		return($arr_degen_score);	
	}
	

sub get_array_degeneracy_score_for_undotted_repeats()
	{
		my($potential_repeat,$ref_arr_repeats)=@_;
		
		my $arr_degen_score=0;
		
		foreach my $repeat(@{$ref_arr_repeats})
			{
				#my $new_r_seq=&change_dots_to_bases($repeat,$model_repeat);
				
				my $sim_score=&get_similarity_score($repeat,$potential_repeat);
				
				my $neg_score=length($potential_repeat)-$sim_score;
				
				$arr_degen_score=$arr_degen_score-$neg_score;
			}
			
		return($arr_degen_score);	
	}


sub search_alternate_repeat_sequence()
	{
		my($range,$accession,$model_repeat,$avg_spacer_length,$current_array,$modified_array)=@_;
		#print "Going to check for consensus sequence with $model_repeat=$model_repeat in $accession:\n\n";
		
		my $case_found=0;
		my $new_model_repeat="";
		my $potential_alternate_repeat="NA";
		
		#---------- now open the sequence file and get the sequence string -------------
		#open(SEQ,"$tmp_dir\/$accession\.fna") or print "$!";
		#my @arr_seq=<SEQ>;
		#close(SEQ);
		#my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;
		
		my @arr_repeats;
		my %hash_of_repeats;
		
		my @arr_model_repeat_bases=split('',$model_repeat);
				
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
				#print "@{$current_array-[0]}\n";
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
							
				#print "C: $current_line\n";
							
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
						
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1];
				
				if(defined $tmp_array[2])
					{
						$s_seq=$tmp_array[2];
					}
				else{
						$s_seq="";
					}
				
				if($tmp_array[3])
					{
						$comment=$tmp_array[3];
					}
				else{
						$comment="";
					}
																				
				$r_length=length($r_seq);
				$s_length=length($s_seq);
				
				
				my $new_r_seq=&change_dots_to_bases($r_seq,$model_repeat);
				
				if($hash_of_repeats{$new_r_seq})
					{
						$hash_of_repeats{$new_r_seq}=$hash_of_repeats{$new_r_seq}+1;
					}
				else{
						$hash_of_repeats{$new_r_seq}=1;
					}
					
					
				push(@arr_repeats,$new_r_seq);		
			}
		
		my @arr_of_model_repeats;
		my $heighest_repeat_occurrence_score=0;
		foreach my $repeat(sort{$hash_of_repeats{$b}<=>$hash_of_repeats{$a}} keys %hash_of_repeats)
			{
				if($hash_of_repeats{$repeat}>1 and $hash_of_repeats{$repeat}>=int(($#arr_repeats+1)*25/100) and $repeat!~/-/ and  $repeat ne $model_repeat)   #check if the alternate repeat covers at least 25% of the total repeats
					{						
						#print "$hash_of_repeats{$repeat}>(($#arr_repeats+1)*25/100)\n";
						$potential_alternate_repeat=$repeat;
						$case_found=1;
						last;
					}	
						
				#$new_model_repeat=$repeat;
				
				#print "\t\t\$new_model_repeat=$new_model_repeat\n";
				#last;
			}
		
		
		#--------- check if there is just one repeat in the array
		#if($#arr_of_model_repeats==0)
		#	{
		#		$new_model_repeat=$arr_of_model_repeats[0];
		#	}
		#else{
		#		$new_model_repeat=$model_repeat;
		#	}
		
		
		
	
			
			#print "\$model_repeat=$model_repeat\n\n";
				
		return($model_repeat,$potential_alternate_repeat,$case_found);		
	}



sub find_array_start_stop_position()
	{
		my ($array_direction,$current_array)=@_;
		
		my $array_start_position=0;;
		my $array_stop_position=0;
		my $average_spacer_length=0;
		my $total_spacer_length=0;	
		my $spacer_count=0;				
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
				#print "@{$current_array-[0]}\n";
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
					
							
				#print "C: $current_line\n";
							
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
						
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1]; my $r_seq_1=$r_seq;$r_seq_1=~s/-//g;
				
				if(not defined $tmp_array[2]){$tmp_array[2]="";}
				$s_seq=$tmp_array[2]; my $s_seq_1=$s_seq;$s_seq_1=~s/-//g;
				
				if($k1<$#{$current_array}-1)
					{
						$total_spacer_length=$total_spacer_length+length($s_seq_1);
						$spacer_count++;
					}
				
				if($tmp_array[3])
					{
						$comment=$tmp_array[3];
					}
				#--------------------------------------------------------
				if($k1==4){$array_start_position=$r_start;}
				if($k1==($#{$current_array}-1))
					{
						if($array_direction !~ "R")
							{
								$array_stop_position=$r_start+length($r_seq_1);
							}
						else{
								$array_stop_position=$r_start-length($r_seq_1);
							}		
					}		
					
			}
		
		if($total_spacer_length>0 and $spacer_count>0)
			{
				$average_spacer_length=int($total_spacer_length/$spacer_count);
			}	
			
		return($array_start_position,$array_stop_position,$average_spacer_length);	
	}


sub get_aligned_region_start_stop()
	{
		
		my($range,$accession,$model_repeat,$ref_repeat)=@_;
		
		my $time1=&get_unique_id();	
				#my $time=$time.int(rand(1000000000)).$arr_letters[rand(int($#arr_letters))].$arr_letters[rand(int($#arr_letters))].$arr_letters[rand(int($#arr_letters))];	

		my $infile=$accession.$range."_".$time1."_model_repeat_bases.tmp";	
		my $outfile=$accession.$range."_".$time1."_clustalw_output.txt";	

		open(WR,">$tmp_dir\/$infile") or print "$!";
		print WR ">M_REPEAT\n$model_repeat\n>R_REPEAT\n$ref_repeat";
		close(WR);
		
						
		my $ret_msg1=`clustalw -QUIET -INFILE=$tmp_dir\/$infile -OUTFILE=$tmp_dir\/$outfile >/dev/null 2>&1`;
		

		#----- read the contents of the output files and store them in arrays, then delete all the files: needed like this, or else files remains undeleted
		open(RD1,"$tmp_dir\/$outfile") or print "$! [$outfile] <br>";
		flock(RD1,2);
		my @arr_rd1=<RD1>;
		close(RD1);


		unlink("$tmp_dir\/$infile") or print "$!\n";
		unlink<$tmp_dir\/*.dnd>;
		unlink("$tmp_dir\/$outfile");	
		


		my $m_repeat_line;
		my $r_repeat_line;
		foreach my $line(@arr_rd1)
			{
				#print $line,"\n";
				chomp $line;$line=~s/\r//g;
				
				if($line=~/^M_REPEAT/)
					{
						my $t_m_repeat_line=$line;
						$t_m_repeat_line=~s/^M_REPEAT//;
						$t_m_repeat_line=~s/^\s+//;
						$m_repeat_line=$m_repeat_line.$t_m_repeat_line
						
					}
				if($line=~/^R_REPEAT/)
					{
						my $t_r_repeat_line=$line;
						$t_r_repeat_line=~s/^R_REPEAT//;
						$t_r_repeat_line=~s/^\s+//;
						$r_repeat_line=$r_repeat_line.$t_r_repeat_line
					}		

			}
		
		return($m_repeat_line,$r_repeat_line);
		
	}
	
sub get_matching_reference_repeat_and_direction()
	{
		
		my($range,$blast_db_file_of_known_repeats,$accession,$model_repeat,$allowed_no_of_mismatches,$lib_of_repeats_with_confirmed_direction)=@_;
		
		#print "\nInside sub get-matching_reference_repeat_and_direction with $accession,$model_repeat :\n";
		
		
		
		
		
		
		my $matching_reference_repeat="NA";
		my $array_direction="NA";
		my $observed_percent_similarity=0;
		
		my $case_found=0;
		my $match_found=0;
		my $ref_repeat_family="NA";
	
		my $allowed_percent_similarity=100-int(($allowed_no_of_mismatches*100)/length($model_repeat)); 				#----- here it has nothing to do with 95% similarity
		my $minimum_length_distribution=100-int(($allowed_no_of_mismatches*100)/length($model_repeat));
		
		

		
		my $u_id=&get_unique_id();
				
				if(length($model_repeat)>=18)
					{
						my $blast_input=$accession.$range.&get_unique_id()."_blast_input.txt";  	open(BIN,">$tmp_dir\/$blast_input");close(BIN);
						my $blast_output=$accession.$range.&get_unique_id()."_blast_output.txt";	open(BOUT,">$tmp_dir\/$blast_output");close(BOUT);
						
						unless(-e "$tmp_dir\/$blast_input")
							{
								#sleep(1);
								select(undef, undef, undef, 0.15);
								open(BIN,">$tmp_dir\/$blast_input");close(BIN);
							}
						
						open(BIN,">$tmp_dir\/$blast_input") or print "$!\n";
						flock(BIN,2);
						print BIN ">MR_1$u_id\n";
						print BIN "$model_repeat";
						close(BIN);
						
						if(-e "$tmp_dir\/$blast_input"){system("chmod 777 $tmp_dir\/$blast_input");}
						#print "Program is going to BLAST the $model_repeat to get matching DRs...\n-db $blast_db_file_of_known_repeats -query $tmp_dir\/$blast_input -out $tmp_dir\/$blast_output \n";				
						system("blastn -task blastn-short -db $blast_db_file_of_known_repeats -query $tmp_dir\/$blast_input -out $tmp_dir\/$blast_output -num_threads 6 -outfmt \"6 sseqid stitle sstart send qseqid qstart qend sstrand length score mismatch gaps sseq qseq\" -word_size 4 >/dev/null 2>&1");							
						#print "\nDone.\n";exit;
						
						if(-e "$tmp_dir\/$blast_output")
							{
								open(BOUT,"$tmp_dir\/$blast_output") or print "$!\n";
								flock(BOUT,2);
								my @arr_bout=<BOUT>;
								close(BOUT);
								
								foreach my $bout_line(@arr_bout)
									{
										chomp $bout_line;$bout_line=~s/\r//g;
										#print "$bout_line\n";
										
										my @arr_line=split('\t',$bout_line); 
										if($#arr_line<12){next;}
										
										my $matching_region_length=$arr_line[8];
										my $mismatches=$arr_line[10];
										my $gaps=$arr_line[11];
												
										my $tmp_matching_ref_repeat_seq=`grep '$arr_line[0]' -A 1 $blast_db_file_of_known_repeats | grep -v '>' >&1`;	
										chomp $tmp_matching_ref_repeat_seq;$tmp_matching_ref_repeat_seq=~s/\r//g;
										
										#print "\$tmp_matching_ref_repeat_seq=$tmp_matching_ref_repeat_seq\n";
												
										if(($matching_region_length-$mismatches-$gaps) >=length($tmp_matching_ref_repeat_seq)*0.90 or ($matching_region_length-$mismatches-$gaps) >=length($model_repeat)*0.90  or ( ($matching_region_length-$mismatches-$gaps) >=length($model_repeat)*0.80) and ($mismatches+$gaps)<=$allowed_no_of_mismatches ) #  (matching_region_length - mismatch)> 90% of MR length or mismatch <2 and length of match >80%
											{
												
												#------ ge the direction --------------------
												if($arr_line[7]=~/minus/)
													{
														$array_direction="R";
													}
												elsif($arr_line[7]=~/plus/)
													{
														$array_direction="F";
													}			
												else{
														$array_direction="Unconfirmed";
													}
												
												#----- get the family ------------------------												
												my @arr_t1=split('_',$arr_line[0]);												
												$ref_repeat_family=$arr_t1[$#arr_t1];
												
												#------ get the original known repeat --------
												my @arr_t2=`grep '$arr_line[0]' -A 1 $blast_db_file_of_known_repeats >&1`;												
												$matching_reference_repeat=$arr_t2[1];
												chomp $matching_reference_repeat;$matching_reference_repeat=~s/\r//g;
												
												#------ get observed_percent_similarity ----------
												
												
												my $length_difference=abs(length($model_repeat)-length($matching_reference_repeat));
												
												if(($mismatches+$gaps+$length_difference) > 0)
													{
														#$observed_percent_similarity=100 - int((($mismatches+$gaps+$length_difference) / length($model_repeat) )*100);
														$observed_percent_similarity=100 - int((($mismatches+$gaps) / $matching_region_length )*100);      ##----- added by ambarish on 11th-Aug-2014
													
													}
												else{
														$observed_percent_similarity=100;
													}		
												
												#print "$bout_line\nMatched: $matching_reference_repeat from family: $ref_repeat_family with direction: $array_direction  and \$observed_percent_similarity=$observed_percent_similarity\n\n";
												
												#$known_repeat_score=3;
												$case_found=1;
												$match_found=1;
												
												last;
											}
									}
								unlink("$tmp_dir\/$blast_output");
							}
						else{
								#print "Error: no BLAST output file found\n";
							}	
						unlink("$tmp_dir\/$blast_input");
					}
				else{
						#$known_repeat_score=-3;
					}
		
		#if($case_found==1)
		#	{
				return($matching_reference_repeat,$ref_repeat_family,$array_direction,$observed_percent_similarity);
		#	}


		#----- following block no longer in use :   by Ambarish on 10/04/2015
		
		#####---- the following block is essential for web-server, as a new DR can be provided by user. without this block, the repeats will never be searched. [added 8th-Aug-2014: ambarish]
		#####
		#------------ first check in the lib of repeats with the model_repeat 
		foreach my $ref_repeat(keys %{$lib_of_repeats_with_confirmed_direction})
			{
				$ref_repeat=~tr/U/T/;

				
						my $model_repeat_rc=$model_repeat; $model_repeat_rc=reverse $model_repeat_rc; $model_repeat_rc=~tr/ACGT/TGCA/;
						
						
						
						if($ref_repeat=~/$model_repeat/ or $model_repeat=~/$ref_repeat/)
							{
								if($ref_repeat=~/$model_repeat/)
									{
										$observed_percent_similarity=int((length($model_repeat)/int(length($ref_repeat)))*100);
									}
								elsif($model_repeat=~/$ref_repeat/)
									{
										$observed_percent_similarity=int((length($ref_repeat)/int(length($model_repeat)))*100);
									}	
								
								if($observed_percent_similarity<80){next;}
									
								$array_direction="F";
								$ref_repeat_family=$lib_of_repeats_with_confirmed_direction->{$ref_repeat};
								
								$match_found=1;
								$case_found=1;
								$matching_reference_repeat=$ref_repeat;
								#print "\n\nForward: [$accession] The $model_repeat matched with $ref_repeat and belongs to group: $lib_of_repeats_with_confirmed_direction->{$ref_repeat} .\n";
								
								#$observed_percent_similarity=100;
								last;
							}
						elsif($ref_repeat=~/$model_repeat_rc/ or $model_repeat_rc=~/$ref_repeat/)
							{
								
								if($ref_repeat=~/$model_repeat_rc/)
									{
										$observed_percent_similarity=int((length($model_repeat_rc)/int(length($ref_repeat)))*100);
									}
								elsif($model_repeat_rc=~/$ref_repeat/)
									{
										$observed_percent_similarity=int((length($ref_repeat)/int(length($model_repeat_rc)))*100);
									}	
								
								if($observed_percent_similarity<80 or $observed_percent_similarity>100){next;}
								
								$array_direction="R";
								$ref_repeat_family=$lib_of_repeats_with_confirmed_direction->{$ref_repeat};
								
								$match_found=1;
								$case_found=1;
								$matching_reference_repeat=$ref_repeat;
								
								#print "\n\nReverse: [$accession] The $model_repeat_rc matched with $ref_repeat and belongs to group: $lib_of_repeats_with_confirmed_direction->{$ref_repeat} .\n";	
								
								#$observed_percent_similarity=100;					
								last;
							}	
							
			}

		if($match_found==0)
			{
				#my $tmp_letters="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
				#my @arr_letters=split('',$tmp_letters);
				
				foreach my $ref_repeat(sort{length($lib_of_repeats_with_confirmed_direction->{$b})<=>length($lib_of_repeats_with_confirmed_direction->{$a})}keys %{$lib_of_repeats_with_confirmed_direction})
					{
						$ref_repeat=~tr/U/T/;
						
						my $top_line="";
						my $bottom_line="";	
				
				
						my $model_repeat_rc=$model_repeat; $model_repeat_rc=reverse $model_repeat_rc; $model_repeat_rc=~tr/ACGT/TGCA/;	
									
						
									
						#while(not $best_bottom_line)
						#{
						
						#----- read the contents of the output files and store them in arrays, then delete all the files: needed like this, or else files remains undeleted
						my @arr_rd1;
						my @arr_rd2;
						
						my $outfile1=&run_water($range,$accession,$model_repeat,$ref_repeat,10,9);	
						
						if(-e "$tmp_dir\/$outfile1")
							{
								open(RD1,"$tmp_dir\/$outfile1") or print "$! [$outfile1] <br>";
								flock(RD1,2);
								@arr_rd1=<RD1>;
								close(RD1);
								
								unlink("$tmp_dir\/$outfile1");
							}
						
						my $outfile2=&run_water($range,$accession,$model_repeat_rc,$ref_repeat,10,9);	
						if(-e "$tmp_dir\/$outfile2")
							{	
								open(RD2,"$tmp_dir\/$outfile2") or print "$! [$outfile2]<br>";
								flock(RD2,2);
								@arr_rd2=<RD2>;
								close(RD2);
								
								unlink("$tmp_dir\/$outfile2") or print "$!\n";
							}
												

							
						
						#------------------- check forward orientation ---------------------------------------------------------------				
						
						
						my $tl_count1=0;
						my $bl_count1=0;
						
						foreach my $line(@arr_rd1)
							{
												#print $line,"\n";
								chomp $line;$line=~s/\r//g;
								if(not $line or $line=~/#/){next;}
								elsif($line=~/\d+/)
									{
										#print "$line\n";
										$line=~s/^\s+//;$line=~s/\s+/\t/g;
										my($start,$seq,$stop)=split('\t',$line);
										if($tl_count1==0){$top_line=$line;$tl_count1++;}
										elsif($bl_count1==0){$bottom_line=$line;$bl_count1++;}
									}						
							}
										
																
						#print "$top_line\n$bottom_line\n\n";	
						
						if($top_line eq "" or $bottom_line eq ""){next;}
						
						my($top_start1,$top_seq1,$top_stop1)=split('\t',$top_line);
						#	$top_start1=$top_start1-1;
						#	$top_stop1=$top_stop1-1;
							
						my($bottom_start1,$bottom_seq1,$bottom_stop1)=split('\t',$bottom_line);
						#	$bottom_start1=$bottom_start1-1;
						#	$bottom_stop1=$bottom_stop1-1;
						
						if(not defined $top_seq1 or $top_seq1=~/-/ or not defined $bottom_seq1 or $bottom_seq1=~/-/){next;}
						
						my $similarity_score1=&get_similarity_score($top_seq1,$bottom_seq1);
						
						#$observed_percent_similarity=$similarity_score1;
						if(length($ref_repeat)>length($model_repeat))
							{
								$observed_percent_similarity=int(($similarity_score1/int(length($ref_repeat)))*100);
							}
						else{
								$observed_percent_similarity=int(($similarity_score1/int(length($model_repeat)))*100);
							}
						#$observed_percent_similarity=int(($similarity_score1/int(length($model_repeat)))*100);
						
						
						
						#print "$ref_repeat ($observed_percent_similarity): $similarity_score1>length($model_repeat)*.95",$similarity_score1>length($model_repeat)*.95,"\n";#exit;
						#if((length($top_seq1)>=length($model_repeat)*0.95 or length($model_repeat)>=length($top_seq1)) and (length($ref_repeat)-$similarity_score1)<=$allowed_no_of_mismatches and $observed_percent_similarity>=80)#>=int(length($ref_repeat)*.95))
						if($observed_percent_similarity<=100 and (length($top_seq1)>=length($model_repeat)*($minimum_length_distribution/100) or length($model_repeat)>=length($top_seq1)) and (length($ref_repeat)-$similarity_score1)<=$allowed_no_of_mismatches and $observed_percent_similarity>=$allowed_percent_similarity)#>=int(length($ref_repeat)*.95))
							{
								$array_direction="F";
								$match_found=1;	
								$case_found=1;
								$matching_reference_repeat=$ref_repeat;
								$ref_repeat_family=$lib_of_repeats_with_confirmed_direction->{$ref_repeat};
								
								
								
								#$observed_percent_similarity=$similarity_score1;
								#$observed_percent_similarity=int(($similarity_score1/int(length($model_repeat)))*100);
								
								#print "\n\nQ: Forward: [$accession] The $model_repeat matched with $ref_repeat and belongs to group: $lib_of_repeats_with_confirmed_direction->{$ref_repeat} .\n";
								#print "\$observed_percent_similarity=$observed_percent_similarity\n";
								#print "\t$top_line\n\t$bottom_line\n\n";
								
								
							}
						else{
								#unlink("$tmp_dir\/$outfile1");
							}	
						
						if($match_found==1)
							{														
								last;
							}
						
						
						
						
						#-------------- now check reverse orientation ---------------------------------------------
						
									
						#---- now open the output file and get the alignment, store the alignment in a compiled file
						
						my $top_line2="";
						my $bottom_line2="";
						
						
						
						my $tl_count2=0;
						my $bl_count2=0;
						
						foreach my $line(@arr_rd2)
							{
												#print $line,"\n";
								chomp $line;$line=~s/\r//g;
								if(not $line or $line=~/#/){next;}
								elsif($line=~/\d+/)
									{
										#print "$line\n";
										$line=~s/^\s+//;$line=~s/\s+/\t/g;
										my($start,$seq,$stop)=split('\t',$line);
										if($tl_count2==0){$top_line2=$line;$tl_count2++;}
										elsif($bl_count2==0){$bottom_line2=$line;$bl_count2++;}
									}						
							}
										
										#print "\n";						
													
						
										
						#print "$top_line\n$bottom_line\n\n";	
						if($top_line2 eq "" or $bottom_line2 eq ""){next;}
						
						my($top_start2,$top_seq2,$top_stop2)=split('\t',$top_line2);
						#	$top_start1=$top_start1-1;
						#	$top_stop1=$top_stop1-1;
							
						my($bottom_start2,$bottom_seq2,$bottom_stop2)=split('\t',$bottom_line2);
						#	$bottom_start1=$bottom_start1-1;
						#	$bottom_stop1=$bottom_stop1-1;
						
						
						if($top_seq2=~/-/ or $bottom_seq2=~/-/){next;}
						if(not defined $top_seq2 or $top_seq2=~/-/ or not defined $bottom_seq2 or $bottom_seq2=~/-/){next;}
						
						my $similarity_score2=&get_similarity_score($top_seq2,$bottom_seq2);
						#$observed_percent_similarity=$similarity_score2;
						if(length($ref_repeat)>length($model_repeat))
							{
								$observed_percent_similarity=int(($similarity_score2/int(length($ref_repeat)))*100);
							}
						else{
								$observed_percent_similarity=int(($similarity_score2/int(length($model_repeat)))*100);
							}		
						
						#if((length($top_seq2)>length($model_repeat)*0.95 or length($model_repeat)>=length($top_seq2))and (length($ref_repeat)-$similarity_score2)<=$allowed_no_of_mismatches and $observed_percent_similarity>=80)#$similarity_score2>=int(length($ref_repeat)*.95))
						if($observed_percent_similarity<=100 and (length($top_seq2)>length($model_repeat)*($minimum_length_distribution/100) or length($model_repeat)>=length($top_seq2)) and (length($ref_repeat)-$similarity_score2)<=$allowed_no_of_mismatches and $observed_percent_similarity>=$allowed_percent_similarity)#$similarity_score2>=int(length($ref_repeat)*.95))
							{					
								
								$array_direction="R";
								$match_found=1;	
								$case_found=1;
								$matching_reference_repeat=$ref_repeat;		
								$ref_repeat_family=$lib_of_repeats_with_confirmed_direction->{$ref_repeat};	
										
								
								
								
								#print "\n\nQ: Reverse: [$accession] The $model_repeat matched with $model_repeat_rc and belongs to group: $lib_of_repeats_with_confirmed_direction->{$ref_repeat} .\n";								
								#print "\$observed_percent_similarity=$observed_percent_similarity\n";
								#print "\t$top_line2\n\t$bottom_line2\n\n";					
								
							}
						else{
								#unlink("$tmp_dir\/$outfile2");
							}
						if($match_found==1){last;}
					}
		
			}
			
			
		return($matching_reference_repeat,$ref_repeat_family,$array_direction,$observed_percent_similarity);	
	}


	
sub get_confidence_score()
	{
		my($total_score,$achieved_score)=@_;
		
		my $confidence="NA";
		
		#-- check if the obtained score is >60% of the total_score -----
		if($achieved_score>($total_score*0.66) and $achieved_score>0.5)   #--0.5 is minimum cutoff to assign any prediction HIGH confidence
			{
				$confidence="HIGH";
			}
		elsif($achieved_score>($total_score*0.33))
			{
				$confidence="MEDIUM";
			}
		else{
				$confidence="LOW";
			}		
		
		return($confidence);
	}	
	
	

sub check_array_direction()
	{
		my($range,$blast_db_file_of_known_repeats,$check_motif_in_repeat,$motif,$check_A_and_T_ratio_in_repeat,$check_similarity_with_reference_repeat,$allowed_no_of_mismatches,$check_secondary_structure_of_repeat,$MFE_cutoff,$MFE_minimum_difference,$MFE_exclude_bases,$check_array_degeneracy,$permitted_mutation_per_array,$check_AT_distribution_in_flanks,$AT_distribution_window,$AT_distribution_minimum_percentage_difference,$check_longer_leader,$Motif_match_score,$A_and_T_ratio_score,$Similarity_score,$MFE_score,$Array_degeneracy_score,$AT_distribution_score,$Longer_leader_score,$array_start_position,$array_stop_position,$accession,$model_repeat,$all_gene_positions_folder,$all_gene_positions_file,$lib_of_repeats_with_confirmed_direction,$current_array,$modified_array)=@_;
		#print "\tGoing to check for consensus sequence with \$model_repeat=$model_repeat and ref:$model_repeats_reference_string for $accession:\n\n";
		
		
		


		my $array_direction="NA";
		my $repeat_family="NA";
		my $matching_reference_repeat="NA";
		my $matching_reference_repeat_direction="NA";
		my $match_found=0;
		my $observed_percent_similarity=0;
			
		#($matching_reference_repeat,$repeat_family,$array_direction,$observed_percent_similarity)=split(';',$model_repeats_reference_string);		
		
		
		#system("echo 'MR:$model_repeat' >log.txt");

		
		
		
		
		#-------------------------- some common parameters -----------------------------------------------------			
		my $case_found=0;
		my $new_model_repeat="";
		my $potential_alternate_repeat="NA";
		my $array_direction_MEMO="";
		

		
		
		#---------- now open the sequence file and get the sequence string -------------------------------------------------
		#system("echo '$accession\.fna' >>log1.txt");
		open(SEQ,"$tmp_dir\/$accession\.fna") or print "$!";
		my @arr_seq=<SEQ>;
		close(SEQ);
		my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;
		#if($species_seq eq "")
		#	{
		#		print "\n\n\n\n\n\n\n\n problem with $tmp_dir\/$accession\.fna\n\n\n\n\n\n\n\n\n\n\n";
		#	}
		
		
		#--------------------------------------------- first check the model_repeats for matching to ref repeat, motif and A_and_T count ----------------------------------------------------------
		
		
		#---- check motif -------------------------------------------------
		if($check_motif_in_repeat==1)
			{					
				
				my $model_repeat_rc=$model_repeat; $model_repeat_rc= reverse $model_repeat_rc; $model_repeat_rc=~tr/ACGT/TGCA/;
				
				my $suggested_direction="NA";
										
				if($model_repeat=~/$motif/ and $model_repeat_rc=~/$motif/){$suggested_direction="NA";}
				elsif($model_repeat=~/$motif/ ){$suggested_direction="F";}
				elsif($model_repeat_rc=~/$motif/){$suggested_direction="R";}

				#-----  change all .? to N
				$motif=~s/\.\?/\(N\)/g;
				
				if($suggested_direction!~/NA/)
					{
						$array_direction_MEMO=$array_direction_MEMO." Motif $motif match prediction:         $suggested_direction Score: $Motif_match_score/$Motif_match_score;";
					}
				else{
						$array_direction_MEMO=$array_direction_MEMO." Motif $motif match prediction:         $suggested_direction Score: 0/$Motif_match_score;";
					}		
				
			}
		else{
				#-----  change all .? to N
				$motif=~s/\.\?/\(N\)/g;
				$array_direction_MEMO=$array_direction_MEMO." Motif $motif matched prediction:       NA ;";
			}	


		#---- check A&T -------------------------------------------------
		if($check_A_and_T_ratio_in_repeat==1)
			{				
				#--- replace all U with T if present -----------
				$model_repeat=~tr/U/T/;

				##################### check direction using no. of As and Ts: mainly applicable to Archaeal genomes #######################		
				my $no_of_As_in_repeat=0;
				my $no_of_Ts_in_repeat=0;
				
				if($model_repeat=~/A/i)
					{								
						$no_of_As_in_repeat=$model_repeat=~s/A/A/gi;
					}
				if($model_repeat=~/T/i)
					{		
						$no_of_Ts_in_repeat=$model_repeat=~s/T/T/gi;	
					}
						
				my $A_and_T_suggestion="NA";
											
				my $at_p_d=0;
				if(($no_of_As_in_repeat+$no_of_Ts_in_repeat)>0 and length($model_repeat)>0)
					{
						$at_p_d=(($no_of_As_in_repeat+$no_of_Ts_in_repeat)/length($model_repeat))*100; $at_p_d=sprintf("%.2f",$at_p_d);
					}
				my $gc_p_d=0;
				if(length($model_repeat)>0)
					{	
						$gc_p_d=((length($model_repeat) -$no_of_As_in_repeat-$no_of_Ts_in_repeat)/length($model_repeat))*100; $gc_p_d=sprintf("%.2f",$gc_p_d);
					}
				#---- skip if AT% is higher than GC% in the repeat -----------------
				if($at_p_d>$gc_p_d)
					{
						$array_direction_MEMO=$array_direction_MEMO." A,T distribution in repeat prediction:     NA [Repeat is AT rich:$at_p_d%AT];";
					}
				else{	
						#------------------------------------------------------------------
											
						if($no_of_As_in_repeat>$no_of_Ts_in_repeat)
							{
								$A_and_T_suggestion="F";
							}
						elsif($no_of_Ts_in_repeat>$no_of_As_in_repeat)
							{
								$A_and_T_suggestion="R";
							}							

				
						$array_direction_MEMO=$array_direction_MEMO." A,T distribution in repeat prediction:     $A_and_T_suggestion [$no_of_As_in_repeat,$no_of_Ts_in_repeat] Score: $A_and_T_ratio_score/$A_and_T_ratio_score;";
					}
			}
		else{
						$array_direction_MEMO=$array_direction_MEMO." A,T distribution in repeat prediction:     NA ;";
			}			
		
		
		#------- check reference-repeats ---------------------------------------------
		if($check_similarity_with_reference_repeat==1)
			{
				
				($matching_reference_repeat,$repeat_family,$array_direction,$observed_percent_similarity)=&get_matching_reference_repeat_and_direction($range,$blast_db_file_of_known_repeats,$accession,$model_repeat,$allowed_no_of_mismatches);
				
				#print "$matching_reference_repeat,$repeat_family,$array_direction,$observed_percent_similarity\n";
				#----------------------------append the suggestion in MEMO string -------------------------------------------	
				if($matching_reference_repeat!~/NA/)
					{
						my $looks_ok=0;
						my $model_repeat_rc=$model_repeat; $model_repeat_rc=reverse $model_repeat_rc; $model_repeat_rc=~tr/ACGTU/TGCAA/;
						
						if($observed_percent_similarity<100 and ($matching_reference_repeat!~/$model_repeat/))
							{
								$looks_ok++;
							}
						if($matching_reference_repeat!~/$model_repeat/)
							{
								$looks_ok++;
							}
						if($model_repeat!~/$matching_reference_repeat/)
							{
								$looks_ok++;
							}
						if($matching_reference_repeat!~/$model_repeat_rc/)
							{
								$looks_ok++;
							}
						if($model_repeat_rc!~/$matching_reference_repeat/)
							{
								$looks_ok++;
							}
										
							
						if($looks_ok>0)
							{
								$matching_reference_repeat_direction=$array_direction;
								$array_direction_MEMO=$array_direction_MEMO." Reference repeat match prediction:         $array_direction [matched $matching_reference_repeat with $observed_percent_similarity\% identity] Score: $Similarity_score/$Similarity_score;";
							}	
					}
				else{
						$array_direction_MEMO=$array_direction_MEMO." Reference repeat match prediction:         NA ;";
					}		
				#-------------------------------------------now check for the AT richness/longer leader/degeneracy -------------------------------
		
			}
		else{
						$array_direction_MEMO=$array_direction_MEMO." Reference repeat match prediction:         NA ;";
			}
		
						
		
	
		#------- check MFE ---------------------------------------------
		if($check_secondary_structure_of_repeat==1)
			{		

				my $forward_strand=$model_repeat;
				my $MFE_suggestion_memo="";
					
				if($MFE_exclude_bases>0)
					{
						my $string_N="";
						for(1..$MFE_exclude_bases){$string_N=$string_N."N";}
						
						$forward_strand=~s/^\S{$MFE_exclude_bases}/$string_N/;
						$forward_strand=~s/\S{$MFE_exclude_bases}$/$string_N/;						
					}	
				
				
				my $other_strand=$forward_strand; $other_strand=reverse $other_strand; $other_strand=~tr/ACGT/TGCA/;
				
				my $mfe_1=&run_rnafold($range,$accession,$forward_strand);				
				my $mfe_2=&run_rnafold($range,$accession,$other_strand);
				
				
				my $MFE_suggestion="NA";	
				
				#---- check if meets minimum valid threshhold ---
				my $valid_case=0;
				if( (abs( abs($mfe_1) - abs($mfe_2) )>=$MFE_minimum_difference) and ( abs($mfe_1)>=$MFE_cutoff or abs($mfe_2)>=$MFE_cutoff))
					{
						$valid_case=1;
					}				
				
				if($mfe_1<$mfe_2 and $valid_case==1){$MFE_suggestion="F";}
				elsif($mfe_2<$mfe_1 and $valid_case==1){$MFE_suggestion="R";}
				

								
				
				if($MFE_suggestion!~/NA/)
					{
						$MFE_suggestion_memo=$MFE_suggestion_memo."$MFE_suggestion [$mfe_1,$mfe_2] Score: $MFE_score/$MFE_score";
					}
				else{
						$MFE_suggestion_memo=$MFE_suggestion_memo."$MFE_suggestion [$mfe_1,$mfe_2] Score: 0/$MFE_score";
					}		
				
				
				#}
				
				$MFE_suggestion_memo=~s/\t$//;
				
				$array_direction_MEMO=$array_direction_MEMO." Secondary Structural analysis prediction:  $MFE_suggestion_memo;";
				#------------------------------------------------------------------------------------------
				#system("echo '$accession\t$array_start_position\t$array_stop_position\t$array_direction [$observed_percent_similarity]\t$MFE_suggestion [$mfe_1,$mfe_2]\t$model_repeat' >>all_repeats_MFE_suggestions.txt");
				#system("echo '$accession\t$array_start_position\t$array_stop_position\t$array_direction [$observed_percent_similarity]\t$MFE_suggestion_memo\t$model_repeat' >>all_repeats_MFE_suggestions.txt");
			}
		else{
				$array_direction_MEMO=$array_direction_MEMO." Secondary structural analysis prediction:  NA ;";
			}		
	
	
	
	
	
	
	
		#------------------- next analyze the other predictions using array information -------------------------------------------------------------------
		

		#-------- check degeneracy -------------------------------------			
		my $degeneracy_suggestion="NA";				
		if($check_array_degeneracy==1 and $#{$current_array}>4)
			{		
						
						#------------------------- check degeneracy -------------------------------------------------------------------
						my $old_top_degeneracy=0;   #keep a backup for printing -----------
						my $old_bottom_degeneracy=0;
						my $top_degeneracy=0;
						my $middle_degeneracy=0;
						my $bottom_degeneracy=0;
						my $degeneracy_in_first_repeat=0;
						my $degeneracy_in_last_repeat=0;
						
						my @arr_repeats;
						my @arr_comments;
							
						for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
							{
								#print "@{$current_array-[0]}\n";
								my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
											
								#print "C: $current_line\n";
											
								my @tmp_array=split('\t',$current_line);
								my($r_start,$r_length,$s_length,$r_seq,$s_seq);
								my $comment="";		
								$r_start=$tmp_array[0];
								$r_seq=$tmp_array[1]; my $r_seq_1=$r_seq;$r_seq_1=~s/-//g;
								
								if(not defined $tmp_array[2]){$tmp_array[2]="";}
								$s_seq=$tmp_array[2]; my $s_seq_1=$s_seq;$s_seq_1=~s/-//g; #if(not defined $s_seq_1){print "\nError with $accession \n\n";}
								
								my $no_of_insertions_and_deletions=0;
								my $dotless_r_seq=$r_seq;$dotless_r_seq=~s/\.+//g;
								
								if($tmp_array[3])
									{
										$comment=$tmp_array[3];$comment=~s/^\s+//;
										if($comment!~/^Del/)
											{
												my @tmp_array=split(' ',$comment);
										
												#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
												my $cur_insertion_bases=$tmp_array[0];
												
												$cur_insertion_bases=~s/,//g;
												$no_of_insertions_and_deletions=length($cur_insertion_bases);
												
												#if(defined $tmp_array[2] and $tmp_array[2]=~/Del/){$no_of_insertions_and_deletions++;}
											}
										else{
												#$no_of_insertions_and_deletions=1;
											}	
									}
								#--------------------------------------------------------							
								push(@arr_repeats,$r_seq);	
								push(@arr_comments,$comment);	
								
								#------ new scoring system ------------------------------
										#my $current_repeat=$r_seq; 
										#my $dotless_current_repeat=$current_repeat; $dotless_current_repeat=~s/\.//g; $$dotless_current_repeat=~s/-+/-/g;
								if($k1==4)
									{
										my $next_line=$$current_array[$k1+1]; chomp $next_line; $next_line=~s/\r+//g; $next_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
										my @tmp_array2=split('\t',$next_line);
										my $next_r_seq=$tmp_array2[1];
										#if($r_seq eq $next_r_seq){$degeneracy_in_first_repeat=$no_of_insertions_and_deletions;}
										#else{
										#		$degeneracy_in_first_repeat=length($current_repeat)+$no_of_insertions_and_deletions;
										#	}
										
										#if($r_seq ne $next_r_seq or $no_of_insertions_and_deletions>0)
										if($dotless_r_seq ne "" or $no_of_insertions_and_deletions>0)
											{
												#print "\tCR:$r_seq\n\tNR:$next_r_seq\n";
												$degeneracy_in_first_repeat=1;
											}
										
									}
								elsif($k1==$#{$current_array}-1)
									{
										my $prev_line=$$current_array[$k1-1]; chomp $prev_line; $prev_line=~s/\r+//g; $prev_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
										my @tmp_array2=split('\t',$prev_line);
										my $prev_r_seq=$tmp_array2[1];
										#if($r_seq eq $prev_r_seq){$degeneracy_in_last_repeat=$no_of_insertions_and_deletions;}
										#else{
										#		$degeneracy_in_last_repeat=length($current_repeat)+$no_of_insertions_and_deletions;
										#	}
										
										
										#if($r_seq ne $prev_r_seq or $no_of_insertions_and_deletions>0){$degeneracy_in_last_repeat=1;}
										if($dotless_r_seq ne "" or $no_of_insertions_and_deletions>0){$degeneracy_in_last_repeat=1;}
									}								
							}
							
							
						my $skip_old_scoring=1;
						if($skip_old_scoring==1)
						{	
							
						my $devider;
						
						if($#arr_repeats<=4){$devider=0.40;}else{$devider=0.33;}   # for shorter arrays which may have only 5 repeats this works well
							
						for(my $i=0;$i<=$#arr_repeats;$i++)
							{
								my $current_repeat=$arr_repeats[$i]; 
								my $dotless_current_repeat=$current_repeat;  $dotless_current_repeat=~s/\.//g; $dotless_current_repeat=~s/-+/-/g;
								my $current_comment;
								my $no_of_insertions_and_deletions=0;
								
								if($arr_comments[$i])
									{
										$current_comment=$arr_comments[$i]; $current_comment=~s/^\s+//g;
										
										if($current_comment!~/^Del/)
											{
												my @tmp_array=split(' ',$current_comment);
										
												#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
												my $cur_insertion_bases=$tmp_array[0];
												
												$cur_insertion_bases=~s/,//g;
												$no_of_insertions_and_deletions=length($cur_insertion_bases);
												
												#if(defined $tmp_array[2] and $tmp_array[2]=~/Del/){$no_of_insertions_and_deletions++;}
											}
										else{
												#$no_of_insertions_and_deletions=1; # do not score deletion from comment,as this will be scored in repeat anyway
											}
									}
								
								if($i<= int($#arr_repeats*$devider))
									{
										#---- check if the current repeat is same as the next repeat or not
										my $prev_r_seq=$arr_repeats[$i-1];
										my $next_r_seq=$arr_repeats[$i+1];
										
										if($i==0) 
											{
												if(length($dotless_current_repeat)==0 or $current_repeat eq $next_r_seq){$top_degeneracy=$top_degeneracy+$no_of_insertions_and_deletions;}
												else{
														$top_degeneracy=$top_degeneracy+length($dotless_current_repeat)+$no_of_insertions_and_deletions;
													}
											}
										else{
												if(length($dotless_current_repeat)==0 or $current_repeat eq $next_r_seq  or $current_repeat eq $prev_r_seq){$top_degeneracy=$top_degeneracy+$no_of_insertions_and_deletions;}
												else{
														$top_degeneracy=$top_degeneracy+length($dotless_current_repeat)+$no_of_insertions_and_deletions;
													}
											}
										
										#$top_degeneracy=$top_degeneracy+length($current_repeat);
										
											#----- check if the first repeat has any degeneracy---
											#if($i==0 and length($current_repeat)>0){$degeneracy_in_first_repeat=1;}
									}
								elsif($i>=int(($#arr_repeats*$devider)*2))
									{
										my $prev_r_seq=$arr_repeats[$i-1];
										my $next_r_seq=$arr_repeats[$i+1];
 

										if($i==$#arr_repeats) 
											{
												if(length($dotless_current_repeat)==0 or $current_repeat eq $prev_r_seq){$bottom_degeneracy=$bottom_degeneracy+$no_of_insertions_and_deletions;}
												else{
														#print "CR:$current_repeat\nPR:$prev_r_seq\n\n";
														$bottom_degeneracy=$bottom_degeneracy+length($dotless_current_repeat)+$no_of_insertions_and_deletions;
													}
											}
										else{
												if(length($dotless_current_repeat)==0 or $current_repeat eq $next_r_seq or $current_repeat eq $prev_r_seq){$bottom_degeneracy=$bottom_degeneracy+$no_of_insertions_and_deletions;}
												else{
														#print "CR:$current_repeat\nPR:$prev_r_seq\n\n";
														$bottom_degeneracy=$bottom_degeneracy+length($dotless_current_repeat)+$no_of_insertions_and_deletions;
													}
											}
										#$bottom_degeneracy=$bottom_degeneracy+length($current_repeat);
											#if($i==$#arr_repeats and length($current_repeat)>0){$degeneracy_in_last_repeat=1;}
									}
								else{
										my $next_r_seq=$arr_repeats[$i+1];
										if(length($dotless_current_repeat)==0 or $current_repeat eq $next_r_seq){$middle_degeneracy=$middle_degeneracy+$no_of_insertions_and_deletions;}
										else{
												$middle_degeneracy=$middle_degeneracy+length($dotless_current_repeat)+$no_of_insertions_and_deletions;
											}
										#$middle_degeneracy=$middle_degeneracy+length($dotless_current_repeat)+$no_of_insertions_and_deletions;
										
									}	
										
								#print "$arr_repeats[$i]\n";
							}

						my $total_array_degeneracy=$top_degeneracy+$middle_degeneracy+$bottom_degeneracy;
						
						#---- add one extra point if the top/bottom first repeat shows degeneracy 
						$old_top_degeneracy=$top_degeneracy;   #keep a backup for printing -----------
						$old_bottom_degeneracy=$bottom_degeneracy;
						$top_degeneracy=$top_degeneracy+$degeneracy_in_first_repeat;			
						$bottom_degeneracy=$bottom_degeneracy+$degeneracy_in_last_repeat;
						#------------------------------------------------------------------------
						
						if($total_array_degeneracy<=$permitted_mutation_per_array){$degeneracy_suggestion="NA";}			
						elsif($bottom_degeneracy>$top_degeneracy){$degeneracy_suggestion="F";}
						elsif($top_degeneracy>$bottom_degeneracy){$degeneracy_suggestion="R";}
						elsif($bottom_degeneracy==$top_degeneracy)
							{
								#if($degeneracy_in_first_repeat!=1 and $degeneracy_in_last_repeat==1){$degeneracy_suggestion="F";}
								#elsif($degeneracy_in_first_repeat==1 and $degeneracy_in_last_repeat!=1){$degeneracy_suggestion="R";}
								#else{
										$degeneracy_suggestion="NA";
								#	}	
							}
						else{$degeneracy_suggestion="NA";}
						}
						
						#----- new scoring -------------------
						#if($degeneracy_in_last_repeat>$degeneracy_in_first_repeat){$degeneracy_suggestion="F";}
						#elsif($degeneracy_in_first_repeat>$degeneracy_in_last_repeat){$degeneracy_suggestion="R";}
						#else{$degeneracy_suggestion="NA";}
						#----------------------------append the suggestion in MEMO string -------------------------------------------	
						if($degeneracy_suggestion!~/NA/)
							{
								$array_direction_MEMO=$array_direction_MEMO." Array degeneracy analysis prediction:      $degeneracy_suggestion [$old_top_degeneracy-$old_bottom_degeneracy] Score: $Array_degeneracy_score/$Array_degeneracy_score;";
							}
						else{
								$array_direction_MEMO=$array_direction_MEMO." Array degeneracy analysis prediction:      $degeneracy_suggestion [$old_top_degeneracy-$old_bottom_degeneracy] Score: 0/$Array_degeneracy_score;";
							}		
						#$array_direction_MEMO=$array_direction_MEMO." $degeneracy_suggestion [$degeneracy_in_first_repeat-$degeneracy_in_last_repeat],";						
						#--------------- record it to a file fore directional analysis ----------------------------------------------
						#system("echo '$accession\t$array_start_position\t$array_stop_position\t$array_direction [$observed_percent_similarity]\t$degeneracy_suggestion [$old_top_degeneracy-$old_bottom_degeneracy]\t$degeneracy_in_first_repeat\t$degeneracy_in_last_repeat' >>array_degeneracy_predictions.txt");
			}
		else{
				$array_direction_MEMO=$array_direction_MEMO." Array degeneracy analysis prediction:      NA ;";
			}					

				
					
		#-------- check AT richness ------------------------------------			
		my $at_distribution="";		
			
		if($check_AT_distribution_in_flanks==1)
			{
				my $at_richness_suggestion="NA";					
				
				
				
						
				#------------------------------- Window for AT richness and at_richness_suggestion -----------------------------------						
						
				my $total_percent_at_in_ls=0;
				my $total_percent_at_in_ts=0;
				
				my $leading_seq="";		
				my $trailing_seq="";
				#-------- now get the AT richness in 150base window before and after the array -------------------------------
				if(($array_start_position-1-$AT_distribution_window)>=0 and ($array_start_position-1-$AT_distribution_window)<length($species_seq))
					{
						$leading_seq=substr($species_seq,($array_start_position-1-$AT_distribution_window),$AT_distribution_window);
					}	
				else{
						$leading_seq=substr($species_seq,0,$array_start_position-1);
					}
					
				if(length($species_seq)>=($array_stop_position-1+$AT_distribution_window))
					{	
						if(($array_start_position-1)>0)
							{
								$trailing_seq=substr($species_seq,$array_stop_position-1,$AT_distribution_window);
							}	
					}
				else{
						if(($array_stop_position-1)<length($species_seq))
							{
								$trailing_seq=substr($species_seq,$array_stop_position-1);
							}
					}		
				my $total_a_in_ls=0;
				my $total_t_in_ls=0;
				
				if(defined $leading_seq and $leading_seq=~/\S+/)
					{				
						$total_a_in_ls=$leading_seq=~s/A/A/gi;
						$total_t_in_ls=$leading_seq=~s/T/T/gi;
					}
				my $total_at_in_ls=$total_a_in_ls+$total_t_in_ls;
				$total_percent_at_in_ls=sprintf("%.1f",($total_at_in_ls/$AT_distribution_window)*100);
								
				my $total_a_in_ts=0;
				my $total_t_in_ts=0;
				
				if(defined $trailing_seq and $trailing_seq=~/\S+/)
					{				
						$total_a_in_ts=$trailing_seq=~s/A/A/gi;
						$total_t_in_ts=$trailing_seq=~s/T/T/gi;
					}	
					
				my $total_at_in_ts=$total_a_in_ts+$total_t_in_ts;
				$total_percent_at_in_ts=sprintf("%.1f",($total_at_in_ts/$AT_distribution_window)*100);
								
				my $higher_richness="NA";
				if(abs($total_percent_at_in_ls-$total_percent_at_in_ts)>=$AT_distribution_minimum_percentage_difference)
					{
					
						if($total_percent_at_in_ls>$total_percent_at_in_ts) # minimum 10% higher
							{
								$higher_richness="F";
							}
						elsif($total_percent_at_in_ts>$total_percent_at_in_ls)
							{
								$higher_richness="R";
							}	
						else{
								$higher_richness="NA";
							}	
					}
				if($higher_richness !~/NA/)
					{				
						$at_distribution=$at_distribution."$higher_richness [$total_percent_at_in_ls-$total_percent_at_in_ts]\%AT Score: $AT_distribution_score/$AT_distribution_score";
					}
				else{
						$at_distribution=$at_distribution."$higher_richness [$total_percent_at_in_ls-$total_percent_at_in_ts]\%AT Score: 0/$AT_distribution_score";
					}				
				#	}# end of for
								
						
				if($total_percent_at_in_ls>$total_percent_at_in_ts)
					{
						$at_richness_suggestion="F";
					}
				elsif($total_percent_at_in_ls<$total_percent_at_in_ts)
					{
						$at_richness_suggestion="R";
					}
				else{
						$at_richness_suggestion="NA";
					}								
						#	}
						
				#----------------------------append the suggestion in MEMO string ------------------------	
				#$array_direction_MEMO=$array_direction_MEMO." $at_richness_suggestion [$total_percent_at_in_ls-$total_percent_at_in_ts]\%AT,";
				$array_direction_MEMO=$array_direction_MEMO." AT richness analysis in flanks prediction: $at_distribution;";
				$at_distribution=~s/\t$//;
						#system("echo '$accession\t$array_start_position\t$array_stop_position\t$array_direction [$observed_percent_similarity]\t$at_distribution' >>at_richness_distribution.txt");		
			}#--- end of AT_richness block						
		else{
				$array_direction_MEMO=$array_direction_MEMO." AT richness analysis in flanks prediction: NA ;";
			}					
		
		
		
		#------ check longer leader ------------------------------------				
		if($check_longer_leader==1 and $all_gene_positions_file ne "NA")
			{	
						
				my $longer_leader_suggestion="NA";
				my $longer_l_suggestions="";	
						
				#---------------------------------- get all the gene positions from the all_gene_and_crispr_positions.txt file --------------------
				#my @arr_gene_positions=`grep -w '$accession' $all_gene_positions_folder/$all_gene_positions_file >&1`; 
				my @arr_gene_positions;#=`grep -w 'CDS' $all_gene_positions_folder/$all_gene_positions_file >&1`;
				&get_all_cds_or_crispr_positions($accession,$all_gene_positions_file,'CDS',\@arr_gene_positions);
				#----------------------------------------------------------------------------------------------------------------------------------
						
						
				if($#arr_gene_positions<=0)
					{
						$array_direction_MEMO=$array_direction_MEMO." Longer leader analysis prediction:         NA ;";
						#next;
					}
				else{

						#------------ now get the previous gene stop and next gene start position -----------------------------------------
						my $previous_gene_stop=0;
						my $next_gene_start=length($species_seq)+1;
						foreach my $gene_det_line(@arr_gene_positions)
								{
									my @arr_tmp1=split('\t',$gene_det_line);
											
									if($arr_tmp1[3]>$previous_gene_stop and $arr_tmp1[3]<$array_start_position)
										{
											$previous_gene_stop=$arr_tmp1[3];
										}
											
									if($arr_tmp1[2]<$next_gene_start and $arr_tmp1[2]>$array_stop_position){$next_gene_start=$arr_tmp1[2];}
								}
						my $left_len=  $array_start_position - $previous_gene_stop + 1;
						my $right_len= $next_gene_start - $array_stop_position + 1;

						
						
						
						
							
						#------------------------------- Longer leader suggestion -----------------------------------------------------------						
						
						my $i=1;
						#for(my $i=0;$i<=1;$i=$i+0.2)  # check and record  the % differences starting with 0% (thats minimum 1 base), 20%, 40%, 60%, 80% and 100%
						#{
							$longer_leader_suggestion="NA";
							if($left_len<$right_len)
								{								
									if($right_len>int($left_len*(1+$i)) and ($left_len>75 or $right_len>75)) #---- there should be a length difference of at least 1.5 times or 75 bases [discuss with chris/peter ]
										{
											$longer_leader_suggestion="R";
										}
									else{
											$longer_leader_suggestion="NA";
										}		
								}
							elsif($left_len>$right_len)
								{
									if($left_len>int($right_len*(1+$i)) and ($left_len>75 or $right_len>75)) #---- there should be a length difference of at least 1.5 times or 75 bases [discuss with chris/peter ]
										{
											$longer_leader_suggestion="F";
										}
									else{
											$longer_leader_suggestion="NA";
										}									
								}
							my $current_diff=(1+$i)*100; $current_diff=$current_diff."%";
							#print "$accession $array_direction [$observed_percent_similarity] |--- $previous_gene_stop ---- $array_start_position ---- $array_stop_position ---- $next_gene_start ----| [$left_len - $right_len] \n";	
							
							if($longer_leader_suggestion!~/NA/)
								{
									$longer_l_suggestions=$longer_l_suggestions."$longer_leader_suggestion [$left_len,$right_len] Score: $Longer_leader_score/$Longer_leader_score";
								}
							else{
									$longer_l_suggestions=$longer_l_suggestions."$longer_leader_suggestion [$left_len,$right_len] Score: 0/$Longer_leader_score";
								}		
							#------------------------------
							
							
						#}
						
						$array_direction_MEMO=$array_direction_MEMO." Longer leader analysis prediction:         $longer_l_suggestions;";
						#system("echo '$accession\t$array_start_position\t$array_stop_position\t$array_direction [$observed_percent_similarity]\tDist: [$left_len-$right_len]\t$longer_l_suggestions' >>longer_leader_prediction.txt");
						
					}
			}
		else{
				$array_direction_MEMO=$array_direction_MEMO." Longer leader analysis prediction:         NA ;";
			}		
						
						
			
		#----------------- now process the MEMO field and get the Final direction ----------------------------------------------------------
		# $A_and_T_ratio_score,$Similarity_score,$MFE_score,$Array_degeneracy_score,$AT_distribution_score,$Longer_leader_score,
		my @arr_memo=split(';',$array_direction_MEMO);
		my $score_F=0;
		my $score_R=0;
		my $total_confidence=0;
		foreach my $analysis(@arr_memo)
			{
				if($analysis=~/Motif/i)
					{
						my $c_score=$Motif_match_score;
											
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1]; $direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
							#$direction=~s/\s+$//;#if($direction=~/ \[\S+\]?$/){$direction=~s/ \[\S+\]?$//;}
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}	
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}					
					}
				elsif($analysis=~/A,T distribution/i)
					{
						my $c_score=$A_and_T_ratio_score;
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}	
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}					
					}
				elsif($analysis=~/Reference repeat/i)
					{
						my $c_score=$Similarity_score;
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}	
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}					
					}
				elsif($analysis=~/Secondary/i)
					{
						my $c_score=$MFE_score;
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}						
					}
				
				
				elsif($analysis=~/Array degeneracy/i)
					{
						my $c_score=$Array_degeneracy_score;
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}		
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}				
					}
				elsif($analysis=~/AT richness/i)
					{
						my $c_score=$AT_distribution_score;
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}						
					}
				elsif($analysis=~/Longer leader/i)
					{
						my $c_score=$Longer_leader_score;					
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}	
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}					
					}						
			}
		
		my $number_of_F=$array_direction_MEMO=~s/F/F/g;
		my $number_of_R=$array_direction_MEMO=~s/R/R/g;
		my $final_array_direction="NA";
		my $final_array_direction_MEMO;	
		
		my $confidence="NA";
		if($score_F>$score_R)
			{
				$final_array_direction="F";
				my $score_diff=$score_F-$score_R;
				$confidence=&get_confidence_score($total_confidence,$score_diff);
			}	
		if($score_F<$score_R)
			{
				$final_array_direction="R";
				my $score_diff=$score_R-$score_F;
				$confidence=&get_confidence_score($total_confidence,$score_diff);
			}
		
		#my $difference_in_prediction=abs(abs($score_F)-abs($score_R));		
		#if($difference_in_prediction>0.66){$confidence="HIGH";}
		#elsif($difference_in_prediction>0.33){$confidence="NORMAL";}
		#else{$confidence="POOR";}
		#$array_direction_MEMO=~s/,$//;
		
		
		$array_direction=$final_array_direction;
		$array_direction_MEMO=$array_direction_MEMO."  ; Final direction:         $final_array_direction [$score_F,$score_R   Confidence: $confidence] ;";
		#print "Final_array_direction\t$final_array_direction\n\n";					
			
		
		
		
		#my $l_flank=substr($species_seq,$array_start_position-1-100,100);
		#my $r_flank=substr($species_seq,$array_stop_position-1,100);
		#-----------------------------------------------------------------------------------------------------------------------------------
		#system("echo 'MR:$model_repeat' >>log.txt");	
		return($matching_reference_repeat,$matching_reference_repeat_direction,$model_repeat,$array_direction,$repeat_family,$array_direction_MEMO,$case_found);	
	
	}



sub revers_an_array()
	{
		my($range,$accession,$model_repeat,$current_array,$modified_array)=@_;
		#print "Going to to reverse the array of $accession:\n\n";
		
		my $case_found=0;
		
				
		#--------- now, as the end of the array is found, now reverse the array
	
		for(my $k1=$#{$current_array}-1;$k1>4;$k1--)
			{
				#print "@{$current_array-[0]}\n";
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
				my $previous_line=$$current_array[$k1-1]; chomp $previous_line; $previous_line=~s/\r+//g; $previous_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
				
							
				#print "C: $current_line\n";		
				
				
				#------------------- process the previous line -------------------------
				
				my @tmp_array2=split('\t',$previous_line);
				my($p_r_start,$p_r_length,$p_s_length,$p_r_seq,$p_s_seq);
				my $p_comment="";		
				$p_r_start=$tmp_array2[0];
				$p_r_seq=$tmp_array2[1]; my $p_r_seq_1=$p_r_seq;$p_r_seq_1=~s/-//g;
				if($tmp_array2[2])
					{
						$p_s_seq=$tmp_array2[2]; 
					}
				else{
						$p_s_seq="";
					}		
				my $p_s_seq_1=$p_s_seq;$p_s_seq_1=~s/-//g;
				my $p_no_of_insertions=0;
				
				if($tmp_array2[3])
					{
						$p_comment=$tmp_array2[3];	$p_comment=~s/^\s+//;
						if($p_comment)
							{					
								if($p_comment=~/^Del/)    # --- deletion anywhere ----
									{
										#$array_total_degeneracy=$array_total_degeneracy+1;
										#next;
									}
								else{										
										my @tmp_arr1=split(' ',$p_comment);
												
												#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
										if(defined $tmp_arr1[0] and defined $tmp_arr1[1])
											{		
												my $p_insertion_bases=$tmp_arr1[0];
												my $p_insertion_positions=$tmp_arr1[1];
													$p_insertion_positions=~s/\[//g;
													$p_insertion_positions=~s/\]//g;
																#print "\$cur_insertion_bases=$cur_insertion_bases\n";
												
												my $comp_p_insertion_bases=$p_insertion_bases;
												$comp_p_insertion_bases=~tr/ACGT/TGCA/;   #--- no need to reverse the string
												
																								
												$p_insertion_bases=~s/,//g;
												$p_no_of_insertions=length($p_insertion_bases);		
												
												$p_comment=$comp_p_insertion_bases." [".$p_insertion_positions."]";
											}
									}
							}
						
					}
				
				
				
				#-------------------- process the current line ---------------------			
				my @tmp_array=split('\t',$current_line);
				my($c_r_start,$c_r_length,$c_s_length,$c_r_seq,$c_s_seq);
				my $c_comment="";		
				$c_r_start=$tmp_array[0];
				$c_r_seq=$tmp_array[1]; my $c_r_seq_1=$c_r_seq;$c_r_seq_1=~s/-//g;
				if(not defined $tmp_array[2]){$tmp_array[2]="";}
				$c_s_seq=$tmp_array[2]; my $c_s_seq_1=$c_s_seq;$c_s_seq_1=~s/-//g;
				
				my $c_no_of_insertions=0;
				
				if($tmp_array[3])
					{
						$c_comment=$tmp_array[3];	$c_comment=~s/^\s+//;
						if($c_comment)
							{					
								if($c_comment=~/^Del/)    # --- deletion anywhere ----
									{
										#$array_total_degeneracy=$array_total_degeneracy+1;
										#next;
									}
								else{										
										my @tmp_arr1=split(' ',$c_comment);
												
												#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
										if(defined $tmp_arr1[0] and defined $tmp_arr1[1])
											{		
												my $c_insertion_bases=$tmp_arr1[0];
												my $c_insertion_positions=$tmp_arr1[1];
													$c_insertion_positions=~s/\[//g;
													$c_insertion_positions=~s/\]//g;
																#print "\$cur_insertion_bases=$cur_insertion_bases\n";
												
												my $comp_c_insertion_bases=$c_insertion_bases;
												$comp_c_insertion_bases=~tr/ACGT/TGCA/;   #--- no need to reverse the string
												
																																	
												$c_insertion_bases=~s/,//g;
												$c_no_of_insertions=length($c_insertion_bases);	
												
												
												$c_comment=$comp_c_insertion_bases." [".$c_insertion_positions."]";	
											}
									}
							}
						
					}
														
				my $new_r_start=$c_r_start+length($c_r_seq_1)+$c_no_of_insertions;		
				my $new_r_seq=$c_r_seq; $new_r_seq=reverse $new_r_seq;$new_r_seq=~tr/ACGT/TGCA/;
				my $new_s_seq=$p_s_seq; $new_s_seq=reverse $new_s_seq;$new_s_seq=~tr/ACGT/TGCA/;
			    #------------------------------------------------------------------------
				#push(@arr_repeats,$r_seq);	
				my $new_rec_line="";
				$new_rec_line="$new_r_start\t$new_r_seq\t$new_s_seq\t$c_comment";
				#print "\tD:$new_rec_line\n";
				push(@{$modified_array},$new_rec_line);
				
				
				#--- now create and push the last line ---------------------------------- 
				if(($k1-1)==4)
					{
						my $l_r_start=$p_r_start+length($p_r_seq_1)+$p_no_of_insertions;
						my $l_r_seq=$p_r_seq;  $l_r_seq=reverse $l_r_seq;$l_r_seq=~tr/ACGT/TGCA/;
						my $l_s_seq="|";
						my $l_comment=$p_comment;
						
						my $l_rec_line="$l_r_start\t$l_r_seq\t$l_s_seq\t$l_comment";
						push(@{$modified_array},$l_rec_line);
					}
					
				
				
				#print "$new_rec_line\n";	
			}
		
		
		$model_repeat=reverse $model_repeat; $model_repeat=~tr/ACGT/TGCA/;
		$case_found=1;
		return($model_repeat,$case_found);
	}





sub forward_a_reversed_array()
	{
		my($accession,$model_repeat,$current_array,$modified_array)=@_;
		#print "Going to to reverse the array of $accession:\n\n";
		
		my $case_found=0;
		
				
		#--------- now, as the end of the array is found, now reverse the array
	
		for(my $k1=$#{$current_array}-1;$k1>4;$k1--)
			{
				#print "@{$current_array-[0]}\n";
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
				my $previous_line=$$current_array[$k1-1]; chomp $previous_line; $previous_line=~s/\r+//g; $previous_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
				
							
				#print "C: $current_line\n";
				
				
				
				
				#------------------- process the previous line -------------------------
				
				my @tmp_array2=split('\t',$previous_line);
				my($p_r_start,$p_r_length,$p_s_length,$p_r_seq,$p_s_seq);
				my $p_comment="";		
				$p_r_start=$tmp_array2[0];
				$p_r_seq=$tmp_array2[1]; my $p_r_seq_1=$p_r_seq;$p_r_seq_1=~s/-//g;
				$p_s_seq=$tmp_array2[2]; my $p_s_seq_1=$p_s_seq;$p_s_seq_1=~s/-//g;
				
				my $p_no_of_insertions=0;
				
				if($tmp_array2[3])
					{
						$p_comment=$tmp_array2[3];	$p_comment=~s/^\s+//;
						if($p_comment)
							{					
								if($p_comment=~/^Del/)    # --- deletion anywhere ----
									{
										#$array_total_degeneracy=$array_total_degeneracy+1;
										#next;
									}
								else{										
										my @tmp_arr1=split(' ',$p_comment);
												
												#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
										my $p_insertion_bases=$tmp_arr1[0];
										my $p_insertion_positions=$tmp_arr1[1];
											$p_insertion_positions=~s/\[//g;
											$p_insertion_positions=~s/\]//g;
														#print "\$cur_insertion_bases=$cur_insertion_bases\n";
										
										my $comp_p_insertion_bases=$p_insertion_bases;
										$comp_p_insertion_bases=~tr/ACGT/TGCA/;   #--- no need to reverse the string
										
																						
										$p_insertion_bases=~s/,//g;
										$p_no_of_insertions=length($p_insertion_bases);		
										
										$p_comment=$comp_p_insertion_bases." [".$p_insertion_positions."]";
									}
							}
						
					}
				
				
				
				#-------------------- process the current line ---------------------			
				my @tmp_array=split('\t',$current_line);
				my($c_r_start,$c_r_length,$c_s_length,$c_r_seq,$c_s_seq);
				my $c_comment="";		
				$c_r_start=$tmp_array[0];
				$c_r_seq=$tmp_array[1]; my $c_r_seq_1=$c_r_seq;$c_r_seq_1=~s/-//g;
				$c_s_seq=$tmp_array[2]; my $c_s_seq_1=$c_s_seq;$c_s_seq_1=~s/-//g;
				
				my $c_no_of_insertions=0;
				
				if($tmp_array[3])
					{
						$c_comment=$tmp_array[3];	$c_comment=~s/^\s+//;
						if($c_comment)
							{					
								if($c_comment=~/^Del/)    # --- deletion anywhere ----
									{
										#$array_total_degeneracy=$array_total_degeneracy+1;
										#next;
									}
								else{										
										my @tmp_arr1=split(' ',$c_comment);
												
												#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
										my $c_insertion_bases=$tmp_arr1[0];
										my $c_insertion_positions=$tmp_arr1[1];
											$c_insertion_positions=~s/\[//g;
											$c_insertion_positions=~s/\]//g;
														#print "\$cur_insertion_bases=$cur_insertion_bases\n";
										
										my $comp_c_insertion_bases=$c_insertion_bases;
										$comp_c_insertion_bases=~tr/ACGT/TGCA/;   #--- no need to reverse the string
										
																															
										$c_insertion_bases=~s/,//g;
										$c_no_of_insertions=length($c_insertion_bases);	
										
										
										$c_comment=$comp_c_insertion_bases." [".$c_insertion_positions."]";	
									}
							}
						
					}
														
				my $new_r_start=$c_r_start-length($c_r_seq_1)-$c_no_of_insertions;		
				my $new_r_seq=$c_r_seq; $new_r_seq=reverse $new_r_seq;$new_r_seq=~tr/ACGT/TGCA/;
				my $new_s_seq=$p_s_seq; $new_s_seq=reverse $new_s_seq;$new_s_seq=~tr/ACGT/TGCA/;
			    #------------------------------------------------------------------------
				#push(@arr_repeats,$r_seq);	
				my $new_rec_line="";
				$new_rec_line="$new_r_start\t$new_r_seq\t$new_s_seq\t$c_comment";
				#print "\tD:$new_rec_line\n";
				push(@{$modified_array},$new_rec_line);
				
				
				#--- now create and push the last line ---------------------------------- 
				if(($k1-1)==4)
					{
						my $l_r_start=$p_r_start-length($p_r_seq_1)-$p_no_of_insertions;
						my $l_r_seq=$p_r_seq;  $l_r_seq=reverse $l_r_seq;$l_r_seq=~tr/ACGT/TGCA/;
						my $l_s_seq="|";
						my $l_comment=$p_comment;
						
						my $l_rec_line="$l_r_start\t$l_r_seq\t$l_s_seq\t$l_comment";
						push(@{$modified_array},$l_rec_line);
					}
					
				
				
				#print "$new_rec_line\n";	
			}
		
		
		$model_repeat=reverse $model_repeat; $model_repeat=~tr/ACGT/TGCA/;
		$case_found=1;
		return($model_repeat,$case_found);
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




sub run_rnafold()
	{
		my($range,$accession,$sequence)=@_;
		my $mfe=0;
		
	
		my $time1=&get_unique_id().$accession.$range.$sequence;	
		my $rep_seq_file=	"mr_".$time1.".txt";				
		my $out_file=		"MFE_".$time1.".mfe";
		
		open(WR,">$tmp_dir\/$rep_seq_file");close(WR);	if(-e "$tmp_dir\/$rep_seq_file"){system("chmod 777 $tmp_dir\/$rep_seq_file");}
		open(WR,">$tmp_dir\/$out_file");close(WR);		if(-e "$tmp_dir\/$out_file"){system("chmod 777 $tmp_dir\/$out_file");}
				
		my $id=">Seq_$time1";				
		my $string=$id."\n".$sequence;
				
		system("echo '$string' >$tmp_dir\/$rep_seq_file");		
		
		
		##---- check RNAfold installation ----------------------------------------------------
		#my $RNAfold=`which RNAfold >&1 2>&1`; chomp $RNAfold; $RNAfold=~s/\r//g;
		#----- run RNAfold -------------------------------------------------------------------		
		system("RNAfold --noLP --noPS < $tmp_dir\/$rep_seq_file >$tmp_dir\/$out_file");       ###### remember apache has no access to /usr/local/bin and the programs in that
								#print "\n$RNAfold --noPS < $tmp_dir\/$rep_seq_file >$tmp_dir\/$out_file\n";
		#-------------------------------------------------------------------------------------		
				
				
		open(MFE,"$tmp_dir\/$out_file") or print "$!: $out_file not found\n";
		flock(MFE,2);
		my @arr_mfe=<MFE>;
		close(MFE); 
				

		foreach my $l(@arr_mfe)
			{						
				chomp $l;$l=~s/\r+//g;
				if($l=~/\)$/)
					{
								#print ": $l\n";
						my @arr_l=split(' ',$l);
						$mfe=$arr_l[$#arr_l];
						$mfe=~s/\(//g;
						$mfe=~s/\)//g;
						$mfe=~s/\s+//g;
								#last;
					}
			}
			
		unlink("$tmp_dir\/$rep_seq_file");
		unlink("$tmp_dir\/$out_file");	
					
		return($mfe);
	}	
	


sub shorten_array()
	{		
		
		my($range,$minimum_no_of_repeats,$allowed_percent_similarity,$user_side_to_shorten,$accession,$model_repeat,$avg_spacer_length,$current_array,$modified_array)=@_;
		
		#print "Going to shorten array of $accession: with $minimum_no_of_repeats,$allowed_percent_similarity,$user_side_to_shorten,$accession,$model_repeat,$avg_spacer_length,$current_array\n\n";
		#system("echo 'US:$user_side_to_shorten' >>log1.txt");
		
		my $case_found=0;
		my $left_flank="";
		my $right_flank="";
		#---------- now open the sequence file and get the sequence string -------------
		#open(SEQ,"$tmp_dir\/$accession\.fna") or print "Error 1: can't find $accession.fna $! <br>\n";
		#my @arr_seq=<SEQ>;
		#close(SEQ);
		#my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;		
		#---- logic ------------------------------------------------------------------------
		#- Step 1: 
		
		
		if($user_side_to_shorten=~/^NA/)
			{
				my $stop_1=0;
				my $stop_2=0;
			
				my $existing_repeats=($#{$current_array}-1)-4;
				
				#print "\$existing_repeats=$existing_repeats\n";
				#-------- first check the 5' end (top to bottom)
				
				for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
					{
						
						my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;							
						#print "A:\$current_line=$current_line\n";
									
						my @tmp_array=split('\t',$current_line);
						my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
								
						$r_start=$tmp_array[0];
						$r_seq=$tmp_array[1];
						if(defined $tmp_array[2])
							{
								$s_seq=$tmp_array[2];
							}
						else{
								$s_seq="";
							}
						
						#$s_seq=$right_flank;
						my $existing_insertion_bases="";
						my $existing_insertion_positions="";
						my $no_of_insertions=0;
						my $number_of_gaps=0;
						$number_of_gaps=$r_seq=~s/-/-/g;		
								
						if($tmp_array[3])
							{						
								$comment=$tmp_array[3];
								
								$comment=~s/^\s+//;
								
								if($tmp_array[3]!~/^\s{0,1}Del/)
									{
										($existing_insertion_bases,$existing_insertion_positions)=split(' ',$comment);
										if($existing_insertion_positions)
											{
												$existing_insertion_positions=~s/\[//g;
												$existing_insertion_positions=~s/\]//g;
											}
											
										$existing_insertion_bases=~s/,//g;	
										$no_of_insertions=length($existing_insertion_bases);	
									}
							}
						else{
								$comment="";
							}	
																						
						#$r_length=length($r_seq);
						#$s_length=length($s_seq);
						
										
						my $repeat_string=&change_dots_to_bases($r_seq,$model_repeat);
						my $similarity_score=&get_similarity_score($repeat_string,$model_repeat);
							$similarity_score=$similarity_score-$number_of_gaps-$no_of_insertions**2;
						
						if($similarity_score<(length($r_seq)*$allowed_percent_similarity/100) and $stop_1==0 and $existing_repeats>$minimum_no_of_repeats)
							{
								$$current_array[$k1]="";
								$existing_repeats--;
								$case_found=1;
							}
						else{
								$stop_1=1;
							}										
					}
					
				#-------- then check the 3' end (bottom)
				
				for(my $k1=$#{$current_array}-1;$k1>=4;$k1--)
					{
						
						my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;	
						
						if($current_line eq ""){next;}		# the line content already deleted from the top in the previous step 
										
						#print "B:\$current_line=$current_line\n";
									
						my @tmp_array=split('\t',$current_line);
						my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
								
						$r_start=$tmp_array[0];
						$r_seq=$tmp_array[1]; #if(not defined $tmp_array[1]){print "ERROR: $current_line";}
						if(defined $tmp_array[2])
							{
								$s_seq=$tmp_array[2];
							}
						else{
								$s_seq="";
							}
						
						#$s_seq=$right_flank;
						my $existing_insertion_bases="";
						my $existing_insertion_positions="";
						my $no_of_insertions=0;
						my $number_of_gaps=0;
						   $number_of_gaps=$r_seq=~s/-/-/g;		
								
						if($tmp_array[3])
							{						
								$comment=$tmp_array[3];
								
								$comment=~s/^\s+//;
								
								if($tmp_array[3]!~/^\s{0,1}Del/)
									{
										($existing_insertion_bases,$existing_insertion_positions)=split(' ',$comment);
										if($existing_insertion_positions)
											{
												$existing_insertion_positions=~s/\[//g;
												$existing_insertion_positions=~s/\]//g;
											}
											
										$existing_insertion_bases=~s/,//g;	
										$no_of_insertions=length($existing_insertion_bases);	
									}
							}
						else{
								$comment="";
							}	
																						
						#$r_length=length($r_seq);
						#$s_length=length($s_seq);
						
										
						my $repeat_string=&change_dots_to_bases($r_seq,$model_repeat);
						my $similarity_score=&get_similarity_score($repeat_string,$model_repeat);
							$similarity_score=$similarity_score-$number_of_gaps-$no_of_insertions**2;
						
						if($similarity_score<(length($r_seq)*$allowed_percent_similarity/100) and $stop_2==0 and $existing_repeats>$minimum_no_of_repeats)
							{
								$$current_array[$k1]="";
								$existing_repeats--;
								$case_found=1;
							}
						else{
								$stop_2=1;
							}										
					}
					
			}
		else{
				
				my $existing_repeats=($#{$current_array}-1)-4;
				
				my@arr_user_sides=split(',',$user_side_to_shorten);
				
				foreach my $side_and_number_of_repeats(@arr_user_sides)
					{
						#print "\$side_and_number_of_repeats=$side_and_number_of_repeats\n";
						my ($side,$number_of_repeats)=split('-',$side_and_number_of_repeats);
						
						if($number_of_repeats>=$existing_repeats){$number_of_repeats=$existing_repeats-1;}
						
						if($side=~/TOP/)
							{
								my $count_1=0;
								for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
									{
										
										if($count_1<$number_of_repeats)
											{
												$$current_array[$k1]="";
												$case_found=1;
											}
										else{
												$existing_repeats=$existing_repeats-$number_of_repeats; #-- update number of remaining repeats
												last;
											}
										$count_1++;	
																				
									}
									
								
							}
						if($side=~/BOTTOM/)
							{
								my $count_2=0;
								for(my $k1=$#{$current_array}-1;$k1>=4;$k1--)
									{
										
										if($count_2<$number_of_repeats)
											{
												$$current_array[$k1]="";
												$case_found=1;
											}
										else{
												$existing_repeats=$existing_repeats-$number_of_repeats; #-- update number of remaining repeats
												last;
											}
										$count_2++;											
									}
							}	
					}
				
			}	
		#---- now remove the blank rows from $$current_array, push them in modified_array -----------
		
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
				
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;	
				if($current_line ne "")
					{
						push(@{$modified_array},$current_line);
					}						
			}
			
		#--- process the last record and replace the last spacer with 0,
		my $last_record=pop(@{$modified_array});
		my @arr_last_record=split('\t',$last_record);	
		$arr_last_record[2]=0;
		my $last_record_new=join("\t",@arr_last_record);
		#---- now push it back ----------------------------------------
		push(@{$modified_array},$last_record_new);	
		
		
			
		return($case_found);		
	}


sub get_repeats_identity()
	{
				my ($range,$accession,$arr_repeats)=@_;
				
				
				my $total_number_of_clusters=0;
				
				my $time_1 = &get_unique_id();		

				$time_1="GSI_".$time_1.$accession.$range;	
				
				my $input_repeat_file=		$time_1."_repeats_in.txt";
				my $output_repeat_file=		$time_1."_repeats_out.txt";
				my $output_cluster_file=	$output_repeat_file.".clstr";
							
				
				open(WR,">$tmp_dir\/$input_repeat_file") or print "$!";	
				
				#foreach my $repeat($$arr_repeats)
				for(my $i=0;$i<=$#{$arr_repeats};$i++)
					{						
						print WR ">S_$i\n$$arr_repeats[$i]\n";
					}
				close(WR);
				
				#print "";
				system("cd-hit-est -i $tmp_dir\/$input_repeat_file -o $tmp_dir\/$output_repeat_file -n 3 -c 0.9 >/dev/null 2>&1");
				
				if(-e "$tmp_dir\/$output_cluster_file")
					{
						$total_number_of_clusters=`grep 'Cluster' $tmp_dir\/$output_cluster_file | wc -l >&1`;
						unlink("$tmp_dir\/$output_cluster_file");
					}	
				
				if(not defined $total_number_of_clusters){$total_number_of_clusters=0;}
				
				unlink("$tmp_dir\/$input_repeat_file");
				unlink("$tmp_dir\/$output_repeat_file");
				
				
				return int($total_number_of_clusters);
				
			}		

sub get_spacers_identity()
	{
				my ($range,$accession,$arr_spacers)=@_;
				
				
				my $total_number_of_clusters=0;
				
				my $time_1 = &get_unique_id();		

				$time_1="GSI_".$time_1.$accession.$range;	
				
				my $input_spacer_file=		$time_1."_spacers_in.txt";
				my $output_spacer_file=		$time_1."_spacers_out.txt";
				my $output_cluster_file=	$output_spacer_file.".clstr";
							
				
				open(WR,">$tmp_dir\/$input_spacer_file") or print "$!";	
				
				#foreach my $spacer($$arr_spacers)
				for(my $i=0;$i<=$#{$arr_spacers};$i++)
					{						
						print WR ">S_$i\n$$arr_spacers[$i]\n";
					}
				close(WR);
				
				#print "";
				system("cd-hit-est -i $tmp_dir\/$input_spacer_file -o $tmp_dir\/$output_spacer_file -n 3 -c 0.8 >/dev/null 2>&1");
				
				if(-e "$tmp_dir\/$output_cluster_file")
					{
						$total_number_of_clusters=`grep 'Cluster' $tmp_dir\/$output_cluster_file | wc -l >&1`;
						unlink("$tmp_dir\/$output_cluster_file");
					}	
				
				if(not defined $total_number_of_clusters){$total_number_of_clusters=0;}
				
				unlink("$tmp_dir\/$input_spacer_file");
				unlink("$tmp_dir\/$output_spacer_file");
				
				
				return int($total_number_of_clusters);
				
			}	
	

sub find_center_position()
	{
		my($coord1,$coord2)=@_;
		my $center_position=0;
		
		my @arr_t2;
		push(@arr_t2,$coord1);
		push(@arr_t2,$coord2);
		@arr_t2=sort{$a<=>$b} @arr_t2;
										
		$center_position=$arr_t2[0]+abs($arr_t2[1]-$arr_t2[0]);
		
		return($center_position);
		
	}


sub create_gff_file()
	{		
		
		my($gff_file,$array_direction,$crispr_index,$array_start_position,$array_stop_position,$accession,$model_repeat,$current_array,$modified_array,$hash_id_lookup_table)=@_;
		
		my($crispr_start_position,$crispr_stop_position)=&sort_two_numbers($array_start_position,$array_stop_position);
		my $crispr_length=$crispr_stop_position-$crispr_start_position+1;
		
		#---------- now open the sequence file and get the sequence string -------------
		open(SEQ1,"$tmp_dir\/$accession\.fna") or print "Can't find $tmp_dir\/$accession\.fna $! \n";
		my @arr_seq=<SEQ1>;
		close(SEQ1);
		my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;
		#---------------------------------------------------------------------------------------
		
		if($hash_id_lookup_table->{$accession})
			{							
				my $original_acc="";
				my $acc_det=$hash_id_lookup_table->{$accession};
				#print "$acc_det\n";
				if($acc_det=~/\|/)
					{
						my @arr_t2=split('\|',$acc_det);
						if($acc_det=~/^gi/)
							{
								$original_acc=$arr_t2[3];
							}
						else{
								$original_acc=$arr_t2[0];
							}	
					}
				elsif($acc_det=~/^\S{1,}-/)
					{
						my @arr_t2=split('-',$acc_det);
						$original_acc=$arr_t2[0];
					}
				elsif($acc_det=~/\s+/)
					{
						my @arr_t2=split(' ',$acc_det);
						$original_acc=$arr_t2[0];
					}						
				else{
						$original_acc=$acc_det;
					}
						
				if($original_acc ne "")
					{
						$accession=$original_acc;
					}	
				#print "$original_acc\n";		
			}
		#---- strand -----------------------
		my $strand="+";
		if($array_direction=~/R/)
			{				
				$strand="-";
			}	
		#-------- first double check the crispr_index 
		my @arr_gff_lines=`grep '$accession' $gff_file | grep 'repeat_region' >&1`;
		if(defined $arr_gff_lines[0] and $arr_gff_lines[0]=~/\S+/)
			{
				$crispr_index=($#arr_gff_lines+1)+1;
			}
		else{
				$crispr_index=1;
			}	
		open(APP,">>$gff_file") or print "$!";
		flock(APP,2);
		
		my $crispr_stop_position_1 = $crispr_stop_position - 1;
		print APP "$accession\tCRISPRDetect\trepeat_region\t$crispr_start_position\t$crispr_stop_position_1\t$crispr_length\t$strand\t.\tID=CRISPR$crispr_index\_$crispr_start_position\_$crispr_stop_position;Note=$model_repeat;Dbxref=SO:0001459;Ontology_term=CRISPR\n";
			
		my $repeat_index=1;	
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
						
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;							
						#print "A:\$current_line=$current_line\n";
									
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
								
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1];
				if(defined $tmp_array[2])
					{
						$s_seq=$tmp_array[2]; 
					}
				else{
						$s_seq="";
					}	
				my $s_seq1=$s_seq;$s_seq1=~s/-//g;
						
				#$s_seq=$right_flank;
				my $existing_insertion_bases="";
				my $existing_insertion_positions="";
				my $no_of_insertions=0;
				my $number_of_gaps=0;
				$number_of_gaps=$r_seq=~s/-/-/g;		
								
				if($tmp_array[3])
					{						
						$comment=$tmp_array[3];
								
						$comment=~s/^\s+//;
								
						if($tmp_array[3]!~/^\s{0,1}Del/)
							{
								($existing_insertion_bases,$existing_insertion_positions)=split(' ',$comment);
								if($existing_insertion_positions)
									{
										$existing_insertion_positions=~s/\[//g;
										$existing_insertion_positions=~s/\]//g;
									}
											
								$existing_insertion_bases=~s/,//g;	
								$no_of_insertions=length($existing_insertion_bases);	
							}
					}
				else{
						$comment="";
					}	
																						

						
										
				
				$comment=~s/Deletion \S+//g;
				
				
			
				#------- repeat ---------------
				#my $repeat_string=&change_dots_to_bases($r_seq,$model_repeat); my $repeat_string1=$repeat_string;$repeat_string1=~s/-//g;
				my $gapless_r_seq=$r_seq;$gapless_r_seq=~s/-//g;				
				#my $repeat_string=&change_dots_to_bases($r_seq,$model_repeat); 
				
				my $repeat_string="";
				if(($r_start-1)>=0 and ($r_start-1)<length($species_seq))
					{
						$repeat_string=substr($species_seq,$r_start-1,length($gapless_r_seq)+$no_of_insertions);				
					}	
				my $repeat_string1=$repeat_string;$repeat_string1=~s/-//g;
				
				my $repeat_length=length($repeat_string);
				
				my $repeat_start=0;	my $repeat_stop=0;
				($repeat_start,$repeat_stop)=&sort_two_numbers($r_start,$r_start+length($repeat_string1)+$no_of_insertions);
	
				my $repeat_stop_1 = $repeat_stop - 1;
				print APP "$accession\tCRISPRDetect\tdirect_repeat\t$repeat_start\t$repeat_stop_1\t$repeat_length\t$strand\t.\tID=CRISPR$crispr_index\_REPEAT$repeat_index\_$repeat_start\_$repeat_stop;Name=CRISPR$crispr_index\_REPEAT$repeat_index\_$repeat_start\_$repeat_stop;Parent=CRISPR$crispr_index\_$crispr_start_position\_$crispr_stop_position;Note=$repeat_string;Dbxref=SO:0001459;Ontology_term=CRISPR\n";
				
				#------- spacer --------------------------------------------------------------------------------------------
				if($k1<$#{$current_array}-1 and length($s_seq1)>0)
					{					
						my $spacer_length=length($s_seq1);			
						my $spacer_start=$repeat_stop;
						my $spacer_stop=$spacer_start+$spacer_length;
						

						($spacer_start,$spacer_stop)=&sort_two_numbers($spacer_start,$spacer_stop);
						
						my $spacer_stop_1 = $spacer_stop - 1;
						print APP "$accession\tCRISPRDetect\tbinding_site\t$spacer_start\t$spacer_stop_1\t$spacer_length\t$strand\t.\tID=CRISPR$crispr_index\_SPACER$repeat_index\_$spacer_start\_$spacer_stop;Name=CRISPR$crispr_index\_SPACER$repeat_index\_$spacer_start\_$spacer_stop;Parent=CRISPR$crispr_index\_$crispr_start_position\_$crispr_stop_position;Note=$s_seq;Dbxref=SO:0001459;Ontology_term=CRISPR\n";	
					}
				#---------------------------------------------------------------------------------------------------------
					
				$repeat_index++;	
										
			}
					
		close(APP);
			
		return($case_found);		
	}

sub sort_two_numbers()
	{
		my($number1,$number2)=@_;
		
		my @arr_t1;
		push(@arr_t1,$number1);
		push(@arr_t1,$number2);
		
		@arr_t1=sort{$a<=>$b} @arr_t1;
		return($arr_t1[0],$arr_t1[1]);
	}



sub calculate_array_quality_score()
	{
		my($range,$species,$accession,$all_gene_positions_folder,$all_gene_positions_file,$matching_reference_repeat,$model_repeat,$current_array)=@_;
		
		my $array_quality_score=0;
		#print "$species,$accession,$all_gene_positions_folder,$all_gene_positions_file,$matching_reference_repeat,$model_repeat\n"; #exit;
		
		#----- first check the presence of Cas1 and Cas2 genes ------------------------------------------------
		my %hash_of_cas_genes;
		my $cas_score=0;
		my $Cas1_found=0;
		my $Cas2_found=0;
		if($all_gene_positions_file ne "NA")
			{
								
				my @arr_cas_genes;	#	=`grep -w 'CRISPR' $all_gene_positions_folder/$all_gene_positions_file >&1`;
				
				&get_all_cds_or_crispr_positions($accession,$all_gene_positions_file,'CRISPR',\@arr_cas_genes);
				
				if(defined $arr_cas_genes[0])
					{
						
						foreach my $cas_gene_det(@arr_cas_genes)
							{
								chomp $cas_gene_det;$cas_gene_det=~s/\r//g;
								my @arr_t1=split('\t',$cas_gene_det);
								
								my $cas_start=$arr_t1[2];
								my $cas_stop=$arr_t1[3];
								my $cas_gene=$arr_t1[4]; $cas_gene=ucfirst($cas_gene);$cas_gene=~s/,//g;
								my $cas_type=$arr_t1[5];$cas_type=~s/,/\//g;
								
								if($cas_gene=~/NA/){next;}
								
								if($cas_gene=~/cas1/i){$Cas1_found=1;}
								if($cas_gene=~/cas2/i){$Cas2_found=1;}
								
								if($Cas1_found==1 and $Cas2_found==1){last;}		
							}
					}
					
				if($Cas1_found==1 or $Cas2_found==1){$cas_score=1;}
				else{
						$cas_score=0;
					}								
			}

		#------------------------------------------------------------------------------------------------------
		
		
		#--------- check if the model repeat or its reverse_comp match any existing repeat from the known repeats ----------
		my $known_repeat_score=0;
		my $model_repeat_rc=$model_repeat; $model_repeat_rc=reverse($model_repeat_rc);$model_repeat_rc=~tr/ACGTU/TGCAA/;
		
		
		if($matching_reference_repeat!~/NA/)
			{
				$known_repeat_score=3;
			}
		else{		
		
				#my $ret_line1="";my $ret_line2="";
				
				#$ret_line1=`grep '$model_repeat' Ref_lib_files/known_valid_repeats.txt >&1`;
				#$ret_line2=`grep '$model_repeat_rc' Ref_lib_files/known_valid_repeats.txt >&1`;
				
				my $u_id=&get_unique_id();
				
				if(length($model_repeat)>=20)
					{
						my $blast_input=$accession.$range.&get_unique_id()."_blast_input.txt";  	open(BIN,">$tmp_dir\/$blast_input");close(BIN);
						my $blast_output=$accession.$range.&get_unique_id()."_blast_output.txt";	open(BOUT,">$tmp_dir\/$blast_output");close(BOUT);
						
						unless(-e "$tmp_dir\/$blast_input")
							{
								sleep(1);
								open(BIN,">$tmp_dir\/$blast_input");close(BIN);
							}
						
						open(BIN,">$tmp_dir\/$blast_input") or print "$!\n";
						flock(BIN,2);
						print BIN ">MR_1$u_id\n";
						print BIN "$model_repeat";
						close(BIN);
						
						if(-e "$tmp_dir\/$blast_input"){system("chmod 777 $tmp_dir\/$blast_input");}
						#print "Program is going to BLAST the $model_repeat to get matching DRs...\n";				
						system("blastn -task blastn-short -db 'DB_KNOWN_REPEATS/known_valid_repeats.fasta' -query $tmp_dir\/$blast_input -out $tmp_dir\/$blast_output -num_threads 6 -outfmt \"6 sseqid stitle sstart send qseqid qstart qend sstrand length score mismatch gaps sseq qseq\" -word_size 4 >/dev/null 2>&1");							
						#print "\nDone.\n";
						
						if(-e "$tmp_dir\/$blast_output")
							{
								open(BOUT,"$tmp_dir\/$blast_output") or print "$!\n";
								flock(BOUT,2);
								my @arr_bout=<BOUT>;
								close(BOUT);
								
								foreach my $bout_line(@arr_bout)
									{
										chomp $bout_line;$bout_line=~s/\r//g;
										#print "$bout_line\n";
										
										my @arr_line=split('\t',$bout_line); 
										if($#arr_line <13){next;}
										if(($arr_line[8]-$arr_line[10]) >=length($model_repeat)*0.90  or ($arr_line[10]<2 and $arr_line[8] >=length($model_repeat)*0.80) ) #  (matching_region_length - mismatch)> 90% of MR length or mismatch <2 and length of match >80%
											{
												#print "$bout_line\n";
												$known_repeat_score=3;
												
												last;
											}
									}
								unlink("$tmp_dir\/$blast_output");
							}
						else{
								#print "Error: no BLAST output file found\n";
							}	
						unlink("$tmp_dir\/$blast_input");
					}
				else{
						#$known_repeat_score=-3;
					}
		}		
		#if($ret_line1 ne "" or $ret_line2 ne "")
		#	{
		#		#print "$ret_line1\n$ret_line2\n";
		#		$known_repeat_score=3;
		#	}
			
			
			
		#------------------------------------------------------------------------------------------------------
		
		#--- check if model_repeat has ATTGAAA(N) in either end -----------------------------------------------
		my $motif_match_score=0;
		if(($model_repeat=~/ATTGAAA/i and $model_repeat=~/\w{15,}ATTGAAA.?$/i) or ($model_repeat_rc=~/ATTGAAA/i and $model_repeat_rc=~/\w{15,}ATTGAAA.?$/i))
			{
				#print "$model_repeat\n$model_repeat_rc\n";
				$motif_match_score=3;
			}
		#------------------------------------------------------------------------------------------------------
		
		
		
		#---- check if overall array identity is >90% ---------------------------------------------------------
		my @arr_spacers;
		my @arr_repeats;
		my @arr_dotted_repeats;
		
		my $overall_repeats_identity_score=0;
		my $avg_percent_identity=0;
		my $nof_deleted_spacers=0;
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
						
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;							
						#print "A:\$current_line=$current_line\n";
									
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
								
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1];
				if(defined $tmp_array[2])
					{
						$s_seq=$tmp_array[2]; 
					}
				else{
						$s_seq="";
					}	
				my $s_seq1=$s_seq;$s_seq1=~s/-//g;				
				my $new_r_seq=&change_dots_to_bases($r_seq,$model_repeat);
				
				push(@arr_dotted_repeats,$r_seq);
				push(@arr_repeats,$new_r_seq);
				
				if($k1<($#{$current_array}-1))
					{
						if($s_seq1 ne "")
							{
								push(@arr_spacers,$s_seq1);
							}
						else{
								push(@arr_spacers,"-");
								$nof_deleted_spacers++;
							}		
					}	
				
				my $similarity_score=&get_similarity_score($new_r_seq,$model_repeat);	#print "\$similarity_score=$similarity_score\n";	
				#print "\$similarity_score=$similarity_score\n";			
				my $percent_identity;
				
				if($similarity_score>0 and length($model_repeat)>0)
					{
						$percent_identity=$similarity_score/length($model_repeat)*100;
						$percent_identity=sprintf("%.1f",$percent_identity);
					}
				   $avg_percent_identity=$avg_percent_identity+$percent_identity;
			}				   
				   
		$avg_percent_identity=$avg_percent_identity/($#arr_repeats+1);$avg_percent_identity=sprintf("%.1f",$avg_percent_identity);
		$overall_repeats_identity_score=sprintf("%.2f",($avg_percent_identity-80)/20);	
			
			
			
			
			
		#---------- check if atleast 2 good repeats present ----------------------------------------------------------------------	
		my $minimum_2_identical_repeats_score=0;
		
		my $identical_repeats_found=0;		
		foreach my $repeat(@arr_dotted_repeats)
			{
				if($repeat!~/[ACGTU-]/i){$identical_repeats_found++;}
			}		
		
		my $total_number_of_dr_clusters=int(&get_repeats_identity($range,$accession,\@arr_repeats)); 
		
		if($identical_repeats_found<2 and $total_number_of_dr_clusters>1)# and $#arr_dotted_repeats>=2)
			{
				$minimum_2_identical_repeats_score=-1.50;
			}
			
		#if($identical_repeats_found<2)# and $#arr_dotted_repeats>=2)
		#	{
		#		$minimum_2_identical_repeats_score=-1.50;
		#	}	

			
		
		#print "\$avg_percent_identity=$avg_percent_identity\t\$overall_repeats_identity_score=$overall_repeats_identity_score\n";








		#----- score the $model_repeat length from the table of distributions ----------------------------------------------------
		my $model_repeat_distribution_score=0;
		open(RD1,"Ref_lib_files/ref_repeat_distribution.txt");
		my @arr_t1=<RD1>;
		close(RD1);
		my %hash_of_repeat_lengths;
		my $highest_peak=0;
		foreach my $line(@arr_t1)
			{
				chomp $line;$line=~s/\r//g;
				if($line=~/#/ or $line eq ""){next;}
				
				my @arr_t2=split('\t',$line);
				
				if($species=~/bacteria/i)
					{
						$hash_of_repeat_lengths{$arr_t2[0]}=$arr_t2[1];
						if($arr_t2[1]>$highest_peak){$highest_peak=$arr_t2[1];}
					}
				elsif($species=~/archaea/i)
					{
						$hash_of_repeat_lengths{$arr_t2[0]}=$arr_t2[2];
						if($arr_t2[2]>$highest_peak){$highest_peak=$arr_t2[2];}
					}	
			}
		
		if(length($model_repeat)>=23 and length($model_repeat)<=47)
			{
				if(defined $hash_of_repeat_lengths{length($model_repeat)})
					{
						$model_repeat_distribution_score=0.25+ sprintf("%.2f", $hash_of_repeat_lengths{length($model_repeat)}/$highest_peak);
					}
				else{
						$model_repeat_distribution_score=0.25;
					}		
			}
		else{
				if(length($model_repeat)<23)
					{
						$model_repeat_distribution_score=-0.25*(23-length($model_repeat));
					}
				elsif(length($model_repeat)>47)
					{
						$model_repeat_distribution_score=-0.25*(length($model_repeat)-47);
					}	
				#$model_repeat_distribution_score=-0.50;
			}	
		
		if($model_repeat_distribution_score< -1){$model_repeat_distribution_score= -1;}
		if($model_repeat_distribution_score> 1){$model_repeat_distribution_score= 1;}
		
		#print "\$model_repeat_distribution_score=$model_repeat_distribution_score\n";
		
		
		
		
		
		
		#------- average spacer length distribution score ------------------------------
		my $average_spacer_length_distribution_score=0;
		my $total_spacer_length_distribution_score=0;

		open(RD2,"Ref_lib_files/ref_spacer_distribution.txt");
		my @arr_ts=<RD2>;
		close(RD2);
		
		my %hash_of_spacer_lengths;
		my $highest_peak2=0;
		foreach my $line(@arr_ts)
			{
				chomp $line;$line=~s/\r//g;
				if($line=~/#/ or $line eq ""){next;}
				
				my @arr_t2=split('\t',$line);
				
				$hash_of_spacer_lengths{$arr_t2[0]}=$arr_t2[1];
				if($arr_t2[1]>$highest_peak2){$highest_peak2=$arr_t2[1];}
	
			}
		foreach my $spacer(@arr_spacers)
			{
				if(length($spacer)>=28 and length($spacer)<=48)
					{
						if(defined $hash_of_spacer_lengths{length($spacer)})
							{
								$total_spacer_length_distribution_score=$total_spacer_length_distribution_score+0.01+ sprintf("%.2f", $hash_of_spacer_lengths{length($spacer)}/$highest_peak2);
							}
						else{
								$total_spacer_length_distribution_score=$total_spacer_length_distribution_score+0.01;
							}		
					}
				elsif(length($spacer)<28)
					{
						my $tmp_spacer_score=-0.10*(28-length($spacer));
						#if($tmp_spacer_score<-1){$tmp_spacer_score= -1;}
						
						$total_spacer_length_distribution_score=$total_spacer_length_distribution_score+$tmp_spacer_score;
					}
				elsif(length($spacer)>48)
					{
						my $tmp_spacer_score=-0.10*(length($spacer)-48);
						#if($tmp_spacer_score<-1){$tmp_spacer_score= -1;}
						
						$total_spacer_length_distribution_score=$total_spacer_length_distribution_score+$tmp_spacer_score;
					}		
			}
			
		if($total_spacer_length_distribution_score!=0)
			{
				if($#arr_spacers>0)
					{
						$average_spacer_length_distribution_score= sprintf("%.2f",$total_spacer_length_distribution_score/$#arr_spacers);
					}
				else{
						$average_spacer_length_distribution_score=$total_spacer_length_distribution_score;
					}		
			}
		if($average_spacer_length_distribution_score< -3){$average_spacer_length_distribution_score = -3;} #--- keeping and lower limit	
		#print "\$average_spacer_length_distribution_score=$average_spacer_length_distribution_score \$total_spacer_length_distribution_score=$total_spacer_length_distribution_score\n";		
		
		
		
		#----------------------------------------------------------------------------------------------------------------------------------------
		my $spacer_identity_score=0;
		#my %hash_of_aligned_spacers;
		#my($spacers_identity)=&check_spacer_identity($accession,\@arr_spacers,\%hash_of_aligned_spacers);
		#$nof_deleted_spacers
		
		my $total_number_of_clusters=int(&get_spacers_identity($range,,$accession,\@arr_spacers));
		
		if(not defined $total_number_of_clusters or $total_number_of_clusters==0)
			{
				select(undef, undef, undef, 1); #--- will sleep for 1/4 seconds
				$total_number_of_clusters=int(&get_spacers_identity($range,$accession,\@arr_spacers));				
			}
		#print "\$total_number_of_clusters= $total_number_of_clusters\n";
		if(defined $total_number_of_clusters and $total_number_of_clusters>0 and $total_number_of_clusters <= int(($#arr_spacers+1 -$nof_deleted_spacers)/2) ){$spacer_identity_score=-3;}
		else{
				$spacer_identity_score=0.2 * $total_number_of_clusters;
				if($spacer_identity_score>1){$spacer_identity_score=1;}
			}
		#print "\$total_number_of_clusters=$total_number_of_clusters <= (($#arr_spacers+1)/2)\n";
		#-----------------------------------------------------------------------------------------------------------------------------------------
		
		#------------------------ score total no of good repeats --------------------------------------------------------------------------------------
		my $no_of_repeats_score=0;
		
		if($identical_repeats_found>0 and $#arr_repeats-$identical_repeats_found>0)
			{
				#$no_of_repeats_score=sprintf("%.2f",log($identical_repeats_found)-log($#arr_repeats-$identical_repeats_found));
				$no_of_repeats_score=sprintf("%.2f",log($#arr_repeats+1)-log($#arr_repeats-$identical_repeats_found));
			}
		elsif($identical_repeats_found>0)
			{
				$no_of_repeats_score=sprintf("%.2f",log($identical_repeats_found));
			}
		if($no_of_repeats_score< -1){$no_of_repeats_score= -1;}
		if($no_of_repeats_score>= 1){$no_of_repeats_score= 1;}
		
				
		#-----------------------------------------------------------------------------------------------------------------------------------------
		my $score_legend="1: cas, 2: likely_repeat, 3: motif_match, 4: overall_repeat_identity, 5: one_repeat_cluster, 6: exp_repeat_length, 7:exp_spacer_length, 8: spacer_identity, 9: log(total repeats) - log(total mutated repeats),";
		my $score_det="1:$cas_score, 2:$known_repeat_score, 3:$motif_match_score, 4:$overall_repeats_identity_score, 5:$minimum_2_identical_repeats_score, 6:$model_repeat_distribution_score, 7:$average_spacer_length_distribution_score, 8:$spacer_identity_score, 9:$no_of_repeats_score,";
		$array_quality_score=$cas_score+$known_repeat_score+$motif_match_score+$overall_repeats_identity_score+$minimum_2_identical_repeats_score+$model_repeat_distribution_score+$average_spacer_length_distribution_score+$spacer_identity_score+$no_of_repeats_score;
		
		return($array_quality_score,$score_det,$score_legend);
		
	}




sub finalize_the_array()
	{
		my($range,$array_quality_score,$score_det,$score_legend,$all_gene_positions_folder,$all_gene_positions_file,$array_direction,$matching_reference_repeat,$repeat_family,$array_direction_MEMO,$potential_alternate_repeat,$accession,$model_repeat,$avg_spacer_length,$left_flank_length,$right_flank_length,$current_array,$modified_array)=@_;
		#print "Going to check for consensus sequence with $model_repeat=$model_repeat in $accession:\n\n";
		
		my $case_found=0;
		my $new_model_repeat="";
		my $left_flank="|";
		my $right_flank="|";
		
		#---------- now open the sequence file and get the sequence string -------------
		open(SEQ1,"$tmp_dir\/$accession\.fna") or print "Can't find $tmp_dir\/$accession\.fna $! \n";
		my @arr_seq=<SEQ1>;
		close(SEQ1);
		my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;
		
		
		
		
		my @arr_repeats;
		my %hash_of_repeats;
		
		my @arr_model_repeat_bases=split('',$model_repeat);
		
		my $longest_spacer_length=0;
		my $average_spacer_length=0;
		my $average_repeat_length=0;
			
		my $array_start_position;
		my $array_stop_position;
							
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
				#print "@{$current_array-[0]}\n";
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
					
							
				#print "C: $current_line\n";
							
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
						
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1]; my $r_seq_1=$r_seq;$r_seq_1=~s/-//g;
				
				if(not defined $tmp_array[2]){$tmp_array[2]="";}
				$s_seq=$tmp_array[2]; my $s_seq_1=$s_seq;$s_seq_1=~s/-//g;
				
				my $rf_nof_insertions=0;
				if($tmp_array[3])
					{
						$comment=$tmp_array[3];
						if($tmp_array[3] and $tmp_array[3]!~/^Del/)
							{
								my $comment2=$tmp_array[3]; chomp $comment2; $comment2=~s/^\s+//;
								#if($comment=~/^Del/){next;}
								
								my @tmp_arr1=split(' ',$comment2);
														
														#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
								my $insertion_bases=$tmp_arr1[0];
								my $insertion_positions=$tmp_arr1[1];
								#print "\t\$insertion_bases=$insertion_bases and \$insertion_positions=$insertion_positions\n";
								if($insertion_bases ne "" and $insertion_bases!~/Del/)
									{																											
										$insertion_bases=~s/,//g;
										$rf_nof_insertions=length($insertion_bases);	
																	
									}
							
								
							}
						
					}
				#--------------------------------------------------------
				if($k1==4){$array_start_position=$r_start;}
				if($k1==($#{$current_array}-1))
					{
						if($array_direction !~ "R")
							{
								$array_stop_position=$r_start+(length($r_seq_1)+$rf_nof_insertions);
							}
						else{
								$array_stop_position = $array_stop_position - 1;
								$array_stop_position=$r_start-(length($r_seq_1)+$rf_nof_insertions);
							}		
					}
				#--------------------------------------------------------																								
				$r_length=length($r_seq_1);
				$s_length=length($s_seq_1);
				
				if($s_length>$longest_spacer_length){$longest_spacer_length=$s_length;}
				
				$average_spacer_length=$average_spacer_length+$s_length;
				$average_repeat_length=$average_repeat_length+$r_length;
							
				#---- check if the spacer length is less than 80% or not and the Flag Del etion is present or not in comment
				if($comment)
					{					
						if($comment=~/Del/)    # --- deletion anywhere ----
							{
								#$array_total_degeneracy=$array_total_degeneracy+1;
								#next;
							}	
					}				
			
				push(@arr_repeats,$r_seq);					
			}
			
		if($#arr_repeats<=0){return('','','','','',$case_found);}
		
		#----- now calculate the $average_spacer_length ------------------------------------------------------------------
		$average_spacer_length=$average_spacer_length/$#arr_repeats;  	$average_spacer_length=sprintf("%.0f",$average_spacer_length);  #--- no need to add one as the last spacer is zero
		$average_repeat_length=$average_repeat_length/($#arr_repeats+1);$average_repeat_length=sprintf("%.0f",$average_repeat_length);
		
		#--- now format the header and footer labels ---------------------------------------------------------------------
		my $label_position   = "  Position";
		my $label_position_u = "==========";
		
		my $label_repeat_len     = "Repeat";
		my $label_repeat_len_u   = "======";
		
		my $label_percent_id   = "   %id";
		my $label_percent_id_u = "======";
		
		my $label_spacer_len     = "Spacer";
		my $label_spacer_len_u   = "======";
		
		my $label_repeat  ="Repeat_Sequence";
		my $label_repeat_u="===============";
		
		if(length($model_repeat)>15)  #--- the minimum length set to 15 so that the structure remains intact
			{
				for(my $i=0;$i<length($model_repeat);$i++)
					{
						if(length($label_repeat)==length($model_repeat)){last;}
						$label_repeat=$label_repeat." ";
						$label_repeat_u=$label_repeat_u."=";
					}
			}
			
		my $label_spacer  ="Spacer_Sequence";	
		my $label_spacer_u="===============";	
		
		if($longest_spacer_length>15)
			{	
				for(my $i=0;$i<$longest_spacer_length;$i++)
					{
						if(length($label_spacer)==$longest_spacer_length){last;}
						$label_spacer=$label_spacer." ";
						$label_spacer_u=$label_spacer_u."=";
					}
			}
		
		my $label_comment    ="Insertion/Deletion";
		my $label_comment_u  ="==================";
		
		my $header_line  ="$label_position\t$label_repeat_len\t$label_percent_id\t$label_spacer_len\t$label_repeat\t$label_spacer\t$label_comment";
		my $header_line_u="$label_position_u\t$label_repeat_len_u\t$label_percent_id_u\t$label_spacer_len_u\t$label_repeat_u\t$label_spacer_u\t$label_comment_u";
		my $footer_line_u="$label_position_u\t$label_repeat_len_u\t$label_percent_id_u\t$label_spacer_len_u\t$label_repeat_u\t$label_spacer_u\t$label_comment_u";
		

		#check if the direction is reverse, mention it
		if($array_direction =~/R/)
			{
				$$current_array[0]=$$current_array[0]."\t\tArray_Orientation: Reverse";
			}
		elsif($array_direction =~/F/)
			{
				$$current_array[0]=$$current_array[0]."\t\tArray_Orientation: Forward";
			}
		else{
				$$current_array[0]=$$current_array[0]."\t\tArray_Orientation: Unconfirmed";
			}		
		
		$$current_array[2]=$header_line;
		$$current_array[3]=$header_line_u;
		$$current_array[$#{$current_array}]=$footer_line_u;
		
		my $avg_percent_identity=0;
		#------------------ now fill the strings with proper number of gaps -------------------------------------------------------------------------------------
		
		my @arr_spacers;
		my $median_spacer_length=&get_median_spacer_length($current_array);
		
		#print "\n\n\$median_spacer_length=$median_spacer_length\n\n\n";
		
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
				#print "@{$current_array-[0]}\n";
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
							
				#print "C: $current_line\n";
				my $next_line="";
				my $next_repeat;
				#my $length_n_repeat=0;
				if($k1<$#{$current_array}-1)
					{
						$next_line=$$current_array[$k1+1]; chomp $next_line; $next_line=~s/\r+//g; $next_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
						my @tmp_array2=split('\t',$next_line);
						$next_repeat=$tmp_array2[1];
						#$next_repeat=~s/-//g;
						#$length_n_repeat=length($next_repeat);
					}
								
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq,$comment);
						
				$r_start=$tmp_array[0];
				$r_seq=$tmp_array[1];my $r_seq1=$r_seq;$r_seq1=~s/-//g;
				
				my $gapless_r_seq1=$r_seq;$gapless_r_seq1=~s/-//g;
				
				if(defined $tmp_array[2])
					{
						$s_seq=$tmp_array[2];
					}	
				else{
						$s_seq="-";
					}
				#--- check if $s_seq has both ATGC and ----
				if($s_seq=~/[ATGC]/i)
					{
						if($s_seq=~/-+$/)
							{
								$s_seq=~s/-+$//;	
							}
						elsif($s_seq=~/^-+/)
							{
								$s_seq=~s/^-+//;
							}		
					}
				elsif($s_seq=~/-+$/)
					{
						$s_seq=~s/-+$/-/;
					}	
					
					
				my $s_seq1=$s_seq;$s_seq1=~s/-//g;
				my $gapless_s_seq1=$s_seq;$gapless_s_seq1=~s/-//g;

				#----- save the spacers in an array for clustalW alignment later
				if($k1< ($#{$current_array}-1))
					{
						if($s_seq1 ne "")
							{
								push(@arr_spacers,$s_seq1);
							}
						else{
								push(@arr_spacers,"-");
							}		
					}
				
				if(defined $tmp_array[3])
					{
						$tmp_array[3]=~s/^\s+//;
						$comment=$tmp_array[3];
					}
				else{
						$comment="";
					}	
																												
				$r_length=length($gapless_r_seq1);
				$s_length=length($gapless_s_seq1);
				
							
				#---- check if the spacer length is less than 80% or not and the Flag Del etion is present or not in comment
				
				#if($r_start=~/1494904/)
				#	{
				#		print "\$s_length=$s_length ($average_spacer_length*80/100)=",($average_spacer_length*80/100),"\n";;
				#	}
				
				if(($comment=~/Del/) and ($k1<$#{$current_array}-1))
					{
						if($s_length>=($median_spacer_length*80/100))
							{
								my @arr_t2=split(' ',$comment);
								for(my $m=0;$m<=$#arr_t2;$m++)
									{
										if($arr_t2[$m]=~/Del/ and $r_seq!~/-$/ and $s_length>=($median_spacer_length*80/100))
											{
												$arr_t2[$m]="";
												$arr_t2[$m+1]="";
											}
									}
								$comment=join(" ",@arr_t2);
								$comment=~s/\s+$//;
									
							}
					}
				elsif(($comment!~/Del/ and $s_length<($median_spacer_length*80/100) and $k1<$#{$current_array}-2 and $next_repeat=~/^-/) or ($comment!~/Del/ and $s_length<($median_spacer_length*80/100) and $k1<$#{$current_array}-1) or ($comment!~/Del/ and $s_length<($median_spacer_length*80/100 and $r_seq=~/-$/)))
					{	
						
						
										
						#if($comment=~/Del/)    # --- deletion anywhere ----
						#	{
						#		#$array_total_degeneracy=$array_total_degeneracy+1;
						#		#next;
						#	}
						my $position;
						if($array_direction=~/R/)
							{
								if($r_seq=~/-$/)
									{
										$position=$r_start-$r_length;
									}
								else{	
										$position=$r_start-$r_length-$s_length;
									}	
							}
						else{
								if($r_seq=~/-$/)
									{
										$position=$r_start+$r_length;
									}
								else{	
										$position=$r_start+$r_length+$s_length;
									}	
							}		
						$comment=$comment." Deletion [$position]";
						$comment=~s/^\s+//;
						#if($r_start=~/1494904/)
						#	{
						#		print "\$s_length=$s_length ($average_spacer_length*80/100)=",($average_spacer_length*80/100),"\$comment=$comment\n";
						#	}	
						
							
					}
				
				
				
				#if($r_start=~/1494904/)
				#	{
				#		print "\$s_length=$s_length ($average_spacer_length*80/100)=",($average_spacer_length*80/100),"\$comment=$comment\n";
				#	}
							
				my $gapped_r_start=&fill_string_with_gaps($r_start,10,"LEFT");
				if($array_direction =~/R/) {
					$gapped_r_start=&fill_string_with_gaps($r_start-1,10,"LEFT");
				}
				
				my $gapless_r_seq=$r_seq;$gapless_r_seq=~s/-//g;
				my $gapped_rep_len=&fill_string_with_gaps(length($gapless_r_seq),6,"LEFT");
				
				my $new_r_seq=&change_dots_to_bases($r_seq,$model_repeat);
				my $similarity_score=&get_similarity_score($new_r_seq,$model_repeat);	#print "\$similarity_score=$similarity_score\n";				
				my $percent_identity=0;
				if($similarity_score>0 and length($model_repeat)>0)
					{
						$percent_identity=$similarity_score/length($model_repeat)*100;
						$percent_identity=sprintf("%.1f",$percent_identity);
						$avg_percent_identity=$avg_percent_identity+$percent_identity;
					}
				my $gapped_percent_identity=&fill_string_with_gaps($percent_identity,6,"LEFT");
				
				my $gapless_s_seq=$s_seq;$gapless_s_seq=~s/-//g;
				my $gapped_spacer_len;
				if($k1==$#{$current_array}-1)
					{
						$gapped_spacer_len=&fill_string_with_gaps(0,6,"LEFT");
						$s_seq="|";
					}
				else{	
						$gapped_spacer_len=&fill_string_with_gaps(length($gapless_s_seq),6,"LEFT");
					}
				
				my $gapped_r_seq;
					if(length($model_repeat)>15)
						{
							$gapped_r_seq=&fill_string_with_gaps($r_seq,length($model_repeat),"RIGHT");
						}
					else{
							$gapped_r_seq=&fill_string_with_gaps($r_seq,15,"RIGHT");
						}	
							
				my $gapped_s_seq;
					if($longest_spacer_length>=15)
						{
							if($s_seq eq ""){$s_seq="-";}
							$gapped_s_seq=&fill_string_with_gaps($s_seq,$longest_spacer_length,"RIGHT");
						}
					else{
							if($s_seq eq ""){$s_seq="-";}
							$gapped_s_seq=&fill_string_with_gaps($s_seq,15,"RIGHT");
						}	
				
				my $gapped_comment="NA";
				if($comment ne ""){$gapped_comment=$comment;}
				$gapped_comment=&fill_string_with_gaps($gapped_comment,7,"RIGHT");
				
				my $updated_line="$gapped_r_start\t$gapped_rep_len\t$gapped_percent_identity\t$gapped_spacer_len\t$gapped_r_seq\t$gapped_s_seq\t$comment";
				#print "$updated_line\n";
				$$current_array[$k1]="$updated_line";
				
				#push(@arr_repeats,$r_seq);		
			}
			
		
		my $gapped_no_of_repeats=&fill_string_with_gaps(($#arr_repeats+1),10,"LEFT");
		my $gapped_avg_repeat_len=&fill_string_with_gaps(length($model_repeat),6,"LEFT");
		
		$avg_percent_identity=$avg_percent_identity/($#arr_repeats+1);$avg_percent_identity=sprintf("%.1f",$avg_percent_identity);
		my $gapped_avg_percent_identity=&fill_string_with_gaps($avg_percent_identity,6,"LEFT");
		
		my $gapped_avg_spacer_len=&fill_string_with_gaps($average_spacer_length,6,"LEFT");
		
		my $gapped_model_repeat=$model_repeat;
		my $gapped_spacer=&fill_string_with_gaps("NA",$longest_spacer_length,"RIGHT");
		
		
		
		my $gapped_comment="NA";
		my $questionable_array=0;
		
		#if($avg_percent_identity<75 or ($average_spacer_length<=10 and length($model_repeat)>=10)){$gapped_comment="Questionable array";$questionable_array=1;}		
		
		$gapped_comment=&fill_string_with_gaps($gapped_comment,7,"RIGHT");				
		
		my $footer_line_misc="$gapped_no_of_repeats\t$gapped_avg_repeat_len\t$gapped_avg_percent_identity\t$gapped_avg_spacer_len\t$gapped_model_repeat\t$gapped_spacer\t$gapped_comment";	
		push(@{$current_array},$footer_line_misc);		
		push(@{$current_array},"");
				
		
		
		
		
		
		
				
		#------ now print the Flanks---------------------------------------------------------------------------------------------------------------------------
		if($array_direction eq "R")
			{
				if(($array_start_position-1+$left_flank_length)<length($species_seq))  #as array_start_position is the highest in the array
					{
						if(($array_start_position-1)<length($species_seq))
							{
								$left_flank=substr($species_seq,$array_start_position-1,$left_flank_length);
							}	
						#print "A: [$array_start_position -1, $left_flank_length] $left_flank\n";
					}
				else{
						if($array_start_position < length($species_seq))
							{
								$left_flank=substr($species_seq,$array_start_position);
								#print "B: [$array_start_position ]", length($species_seq),"\n";
							}						
					}
					
				if(($array_stop_position-1-$right_flank_length)>=0)
					{
						$right_flank=substr($species_seq,($array_stop_position-1-$right_flank_length),$right_flank_length);
						#print "C: [($array_stop_position-1-$right_flank_length),$right_flank_length)] $right_flank\n";
					}
				else{
						if($array_stop_position>0)
							{
								$right_flank=substr($species_seq,0,$array_stop_position);
								#print "D: [$array_stop_position], $right_flank\n";
							}	
					}
				
				if(defined $left_flank and $left_flank=~/\S+/)
					{
						$left_flank=reverse $left_flank;   $left_flank=~tr/ACGT/TGCA/;				
					}	
				if(defined $right_flank and $right_flank=~/\S+/)
					{	
						$right_flank=reverse $right_flank; $right_flank=~tr/ACGT/TGCA/;
					}	
			}
		else{
				#------ now print the Flanks---------------------------------------------------------------------------------------------------------------------------
				if(($array_start_position-$left_flank_length)>0)
					{
						
						$left_flank=substr($species_seq,($array_start_position-1-$left_flank_length),$left_flank_length);
						#print "A: [$array_start_position - $left_flank_length] $left_flank\n";
					}
				else{
						if($array_start_position>1)
							{
								$left_flank=substr($species_seq,0,($array_start_position));
								#print "B: [$array_start_position - 1] $left_flank\n";
							}						
					}
					
				if(($array_stop_position-1+$right_flank_length)<=length($species_seq))
					{
						$right_flank=substr($species_seq,$array_stop_position-1,$right_flank_length);
						#print "C: $right_flank\n";
					}
				else{
						if($array_stop_position<length($species_seq))
							{
								$right_flank=substr($species_seq,$array_stop_position);								
							}
						#print "D: $right_flank\n",length($species_seq)," \$array_stop_position=$array_stop_position\n";
					}
			}
		#-----------------------------------------------------------------------------------------------------------
		if(not defined $left_flank or $left_flank eq "" or length($left_flank)>20000 ){$left_flank="|";}
		if(not defined $right_flank or $right_flank eq "" or length($right_flank)>20000){$right_flank="|";}
		
		#push(@{$current_array},"");	
		push(@{$current_array},"# Left flank :   $left_flank");
		push(@{$current_array},"# Right flank :  $right_flank");
		#push(@{$current_array},"");
		push(@{$current_array},"");	
		
		
		
		#------------------------------------ check if questionable_array------------------------------------------------------------------
		#------ check if the spacers are >=80% identical with cd- hit-est -------
		
		

			
		my $total_number_of_clusters=&get_spacers_identity($range,$accession,\@arr_spacers);
		
		#-----------------------------------------------------------------------
		if($array_quality_score < 1.50)
			{
				$questionable_array=1;
			}
		
		elsif($total_number_of_clusters==1 and $#arr_spacers>=1)
			{
				$questionable_array=1;
				#push(@{$current_array},"# Questionable array : YES [Potential tandem repeat]\t Score: $array_quality_score");
			}
		elsif(length($model_repeat)>=70 and $average_spacer_length<=int(length($model_repeat)/2))
			{
				$questionable_array=1;
			}	
		
			
		if($questionable_array==1)
			{
				if($total_number_of_clusters==1 and $#arr_spacers>=1)
					{
						push(@{$current_array},"# Questionable array : YES [Potential tandem repeat]\t Score: $array_quality_score");	
					}
				elsif($average_spacer_length<=10 and length($model_repeat)>=10) #----- tandem repeats
					{
						push(@{$current_array},"# Questionable array : YES [Potential tandem repeat]\t Score: $array_quality_score");
					}		
				elsif($avg_percent_identity<75)
					{
						push(@{$current_array},"# Questionable array : YES [Average repeats identity <75%]\t Score: $array_quality_score");
					}
				else{
						push(@{$current_array},"# Questionable array : YES\t Score: $array_quality_score");
					}	
						
			}
				
		else{
				push(@{$current_array},"# Questionable array : NO\t Score: $array_quality_score");
			}
			
			
		push(@{$current_array},"# \tScore Detail : $score_det");	
		push(@{$current_array},"# \tScore Legend : $score_legend");	
				
		#-------------------------------------------------------------------------------------------------------
		
		

		
		push(@{$current_array},"# Primary repeat :     $model_repeat");
		push(@{$current_array},"# Alternate repeat :   $potential_alternate_repeat");
		push(@{$current_array},"");
		
		#---- array direction ---------------------------------------------------------
		#push(@{$current_array},"# Array orientation : $array_direction [Ref: $matching_reference_repeat, Repeat_family: $repeat_family ]");
		
		#push(@{$current_array},"# Orientation analysis summary: [$array_direction_MEMO]");
		
		#print FPA "# Directional analysis summary from each method: \n";
		if($array_direction_MEMO ne "NA")
		{
		push(@{$current_array},"# Directional analysis summary from each method:");
		
				my @arr_memo=split(';',$array_direction_MEMO);
				foreach my $memo(@arr_memo)
					{
						$memo=~s/\t+/, /g;$memo=~s/^\s+//;
						#print FPA "# $memo\n";
						
						if($memo eq ""){$memo="----------------------------------------------------------------------------";}
						push(@{$current_array},"# \t$memo");
					}				
				#print FPA "\n";	
		push(@{$current_array},"");
		}
		
		
		
		#-------- print the CAS genes position(s) -----------------------------------------------------------------------------


			
		
		my %hash_of_cas_genes;
		#print "$all_gene_positions_file\n\n";
		my $identified_cas_genes="";
		if($all_gene_positions_file ne "NA")
			{
				
				my $crispr_center_position=&find_center_position($array_start_position,$array_stop_position);
				
				my @arr_cas_genes;			#=`grep -w 'CRISPR' $all_gene_positions_folder/$all_gene_positions_file >&1`;
				
				&get_all_cds_or_crispr_positions($accession,$all_gene_positions_file,'CRISPR',\@arr_cas_genes);
				
				if(defined $arr_cas_genes[0])
					{
						
						foreach my $cas_gene_det(@arr_cas_genes)
							{
								chomp $cas_gene_det;$cas_gene_det=~s/\r//g;
								my @arr_t1=split('\t',$cas_gene_det);
								
								my $cas_start=$arr_t1[2];
								my $cas_stop=$arr_t1[3];
								
								#------- check if the Cas gene is within 25000 upstream or downstream of the array ---------
								my $current_cas_gene_center=&find_center_position($cas_start,$cas_stop);
								#if(abs($current_cas_gene_center-$crispr_center_position)>50000){next;}
								
								
								my $cas_gene=$arr_t1[4]; $cas_gene=ucfirst($cas_gene);$cas_gene=~s/,//g;
								my $cas_protein_acc=$arr_t1[5];$cas_protein_acc=~s/,/\//g;
								
								if($cas_gene=~/NA/){next;}
								
								my $orig_cas_gene_def=$cas_gene;
								if($cas_gene=~/-/)
									{
										my @arr_t2=split('-',$cas_gene);
										$cas_gene=$arr_t2[0];
									}
								if($hash_of_cas_genes{$cas_gene})
									{
										#---- keep the closest one ---------------
										
										my $existing_cas_gene_det=$hash_of_cas_genes{$cas_gene};
										if($existing_cas_gene_det=~/ \[(\d+)-(\d+)\]/)
											{
												
												my $pre_cas_gene_center=&find_center_position($1,$2);												
														#my $current_cas_gene_center=&find_center_position($cas_start,$cas_stop);
												
												#if(abs($crispr_center_position-$pre_cas_gene_center)>abs($crispr_center_position-$current_cas_gene_center))
												#	{
														$hash_of_cas_genes{$cas_gene}="$orig_cas_gene_def:$cas_protein_acc [$cas_start-$cas_stop], ";
												#	}
											}
										
														#$hash_of_cas_genes{$cas_gene}=$hash_of_cas_genes{$cas_gene}."$cas_gene [$cas_start-$cas_stop], ";
									}
								else{
										$hash_of_cas_genes{$cas_gene}="$orig_cas_gene_def:$cas_protein_acc [$cas_start-$cas_stop], ";
									}	
								
		
							}
					}
				
				
				my %valid_cas_genes;
				$valid_cas_genes{'Cas1'}=1;
				$valid_cas_genes{'Cas2'}=1;
				$valid_cas_genes{'Cas3'}=1;
				$valid_cas_genes{'Cas4'}=1;
				$valid_cas_genes{'Cas5'}=1;
				$valid_cas_genes{'Cas6'}=1;
				$valid_cas_genes{'Cas7'}=1;
				$valid_cas_genes{'Cas8'}=1;
				$valid_cas_genes{'Cas9'}=1;
				$valid_cas_genes{'Cas10'}=1;
				$valid_cas_genes{'Cse1'}=1;
				$valid_cas_genes{'Cse2'}=1;
				$valid_cas_genes{'Cas6e'}=1;
				$valid_cas_genes{'Cmr1'}=1;
				$valid_cas_genes{'Cmr3'}=1;
				$valid_cas_genes{'Cmr4'}=1;
				$valid_cas_genes{'Cmr5'}=1;
				$valid_cas_genes{'Cmr6'}=1;
				
				$valid_cas_genes{'Cas8a1'}=1; $valid_cas_genes{'Csa8a1'}=1;$valid_cas_genes{'Csa8a2'}=1;
				$valid_cas_genes{'Cas8b'}=1;
				$valid_cas_genes{'Cas8c'}=1;
				$valid_cas_genes{'Cas10d'}=1;
				$valid_cas_genes{'Cse1'}=1;
				$valid_cas_genes{'Csy1'}=1;
				
				$valid_cas_genes{'Csm2'}=1;
				$valid_cas_genes{'Csn2'}=1;
				
				#my $identified_cas_genes="";				
				foreach my $cas_gene(sort keys %hash_of_cas_genes)
					{
						$cas_gene=ucfirst($cas_gene);
						#---- only match the first part of the Cas genes by removing the number, so Cas6e will check for only Cas
						my @arr_t1=split('\d',$cas_gene);
						my $matched_vcg=0;
						foreach my $valid_cg(keys %valid_cas_genes)
							{
								if($valid_cg=~/$arr_t1[0]/i){$matched_vcg=1;}
							}
						
						#if($matched_vcg==1)
						#	{								
								$hash_of_cas_genes{$cas_gene}=~s/, $//;
								$identified_cas_genes=$identified_cas_genes."$hash_of_cas_genes{$cas_gene}; ";
						#	}
					}
					
					
				if($identified_cas_genes!~/\S+/)
					{
						$identified_cas_genes="NA";
					}
				push(@{$current_array},"# Identified Cas genes:  $identified_cas_genes");
				
			}
		
		#------- array family ---------------------------------------------------------------------------------------------
		my $potential_repeat_family="";
		
		if($repeat_family!~/NA/)
			{
				#push(@{$current_array},"# Array family :         $repeat_family  [Matched known repeat from this family]");
				$potential_repeat_family=$potential_repeat_family."$repeat_family [Matched known repeat from this family], ";
			}
		#else{
				
		my $crispr_family_found_by_cas_genes=0;		
				#----- check the Cas genes of Type-I---------------------------------------------
				if($hash_of_cas_genes{'Cas3'} and ($hash_of_cas_genes{'Cas8a1'} or $hash_of_cas_genes{'Cas8a2'}) )
					{
						if($hash_of_cas_genes{'Cas8a1'})
							{
								$potential_repeat_family=$potential_repeat_family."I-A [Cas3,Cas8a1]/";
								$crispr_family_found_by_cas_genes++;
							}
						elsif($hash_of_cas_genes{'Cas8a2'})
							{
								$potential_repeat_family=$potential_repeat_family."I-A [Cas3,Cas8a2]/";
								$crispr_family_found_by_cas_genes++;
							}		
					}
				if($hash_of_cas_genes{'Cas3'} and ($hash_of_cas_genes{'Csa8a1'} or $hash_of_cas_genes{'Csa8a2'}) )
					{
						if($hash_of_cas_genes{'Csa8a1'})
							{
								$potential_repeat_family=$potential_repeat_family."I-A [Cas3,Csa8a1]/";
								$crispr_family_found_by_cas_genes++;
							}
						elsif($hash_of_cas_genes{'Csa8a2'})
							{
								$potential_repeat_family=$potential_repeat_family."I-A [Cas3,Csa8a2]/";
								$crispr_family_found_by_cas_genes++;
							}
					}
						
				if($hash_of_cas_genes{'Cas3'} and $hash_of_cas_genes{'Cas8b'})
					{
						$potential_repeat_family=$potential_repeat_family."I-B [Cas3,Cas8b]/";
						$crispr_family_found_by_cas_genes++;
					}
				if($hash_of_cas_genes{'Cas3'} and $hash_of_cas_genes{'Cas8c'})
					{
						$potential_repeat_family=$potential_repeat_family."I-C [Cas3,Cas8c]/";
						$crispr_family_found_by_cas_genes++;
					}
				if($hash_of_cas_genes{'Cas3'} and $hash_of_cas_genes{'Cas10d'})
					{
						$potential_repeat_family=$potential_repeat_family."I-D [Cas3,Cas10d]/";
						$crispr_family_found_by_cas_genes++;
					}
				if($hash_of_cas_genes{'Cas3'} and $hash_of_cas_genes{'Cse1'})
					{
						$potential_repeat_family=$potential_repeat_family."I-E [Cas3,Cse1]/";
						$crispr_family_found_by_cas_genes++;
					}				
				if($hash_of_cas_genes{'Cas3'} and $hash_of_cas_genes{'Csy1'})
					{
						$potential_repeat_family=$potential_repeat_family."I-F [Cas3,Csy1]/";
						$crispr_family_found_by_cas_genes++;
					}
					
				#------- type II ------------------------------------------------	
				if($hash_of_cas_genes{'Cas9'} and $hash_of_cas_genes{'Csn2'})
					{
						$potential_repeat_family=$potential_repeat_family."II-A [Cas9,Csn2]/";
						$crispr_family_found_by_cas_genes++;
					}
				if($hash_of_cas_genes{'Cas9'} and $hash_of_cas_genes{'Cas4'})
					{
						$potential_repeat_family=$potential_repeat_family."II-B [Cas9,Cas4]/";
						$crispr_family_found_by_cas_genes++;
					}

				#------- type III ------------------------------------------------	
				if($hash_of_cas_genes{'Cas10'} and $hash_of_cas_genes{'Csm2'})
					{
						$potential_repeat_family=$potential_repeat_family."III-A [Cas10,Csm2]/";
						$crispr_family_found_by_cas_genes++;
					}
				if($hash_of_cas_genes{'Cas10'} and $hash_of_cas_genes{'Cmr5'})
					{
						$potential_repeat_family=$potential_repeat_family."III-B [Cas10,Cmr5]/";
						$crispr_family_found_by_cas_genes++;
					}

					
				#-------------------------------------------------------------------
				$potential_repeat_family=~s/\/$//;
				if($potential_repeat_family eq ""){$potential_repeat_family="NA";}
				if($potential_repeat_family eq "NA")
					{
						push(@{$current_array},"# Array family : NA");
					}
				else{	
						if(defined $identified_cas_genes and $identified_cas_genes=~/\S+/ and $crispr_family_found_by_cas_genes>0)
							{
								push(@{$current_array},"# Array family : $potential_repeat_family [Prediction based on the presence of specific Cas genes]");	
							}
						else{
								push(@{$current_array},"# Array family : $potential_repeat_family  ");
							}	
					}	
										
		#	}
		
		#----- check if the strain info is provided ---------
		if($all_gene_positions_file ne "NA")
			{
				my $strain="NA";
				my @ret=`grep -w '$accession' $tmp_dir\/$all_gene_positions_file | grep 'STRAIN' >&1`;
				if(defined $ret[0] and $ret[0]=~/\S+/)
					{						
						my @arr_s1=split('\t',$ret[0]);
						$strain=$arr_s1[4]; chomp $strain; $strain=~s/\r//;
						push(@{$current_array},"# Sequence source strain : $strain");
					}
			}
		#----- check if the TAXONOMY info is provided ---------
		if($all_gene_positions_file ne "NA")
			{
				my $taxonomy="NA";
				my @ret=`grep -w '$accession' $tmp_dir\/$all_gene_positions_file | grep 'TAXONOMY' >&1`;
				if(defined $ret[0] and $ret[0]=~/\S+/)
					{						
						my @arr_s1=split('\t',$ret[0]);
						$taxonomy=$arr_s1[5]; chomp $taxonomy; $taxonomy=~s/\r//;
						push(@{$current_array},"# Taxonomy hierarchy : $taxonomy");
					}
			}	
		#----------------------------------------------------
		push(@{$current_array},"//");	
		#--------------------------------------------------------------------------------------------------------------------------------------------------------
		$case_found=1;
			
			#print "\$model_repeat=$model_repeat\n\n";
				
		return($model_repeat,$potential_alternate_repeat,$questionable_array,$array_start_position,$array_stop_position,$case_found);		
	}



sub get_unique_id()
	{
		my $letters="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
		my @arr_letters=split('',$letters);		
		my $unique_id = time();	   
		for(my $i = 0; $i < int(rand(50)); $i++)
			 {
				$unique_id .= $arr_letters[int(rand($#arr_letters))];
			 }
		return $unique_id;	 
	}


sub get_median_spacer_length()
	{
		my $current_array=shift(@_);
		
		my $median_spacer_length=0;
		
		my @arr_spacer_lengths;
		for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
			{
				#print "@{$current_array-[0]}\n";
				my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
					
				
				my @tmp_array=split('\t',$current_line);
				
				if(not defined $tmp_array[2]){$tmp_array[2]="";}
				$s_seq=$tmp_array[2]; my $s_seq_1=$s_seq;$s_seq_1=~s/-//g;
				
				push(@arr_spacer_lengths,length($s_seq_1));
			}		
			
		@arr_spacer_lengths=sort{$a<=>$b} @arr_spacer_lengths;
		
		my $number_of_spacers=$#arr_spacer_lengths+1;
				
		if($number_of_spacers % 2 == 0) #--- even number of spacers
			{
				#print "Even number_of_spacers\n";
				my $middle_pos1=sprintf("%.0f",($number_of_spacers/2));				
				my $middle_pos2=sprintf("%.0f",($number_of_spacers/2))+1;
				
				$median_spacer_length=int(($arr_spacer_lengths[int($middle_pos1-1)] + $arr_spacer_lengths[int($middle_pos2-1)])/2);
				
			}
		else{
				my $middle_pos=sprintf("%.0f",($number_of_spacers/2));
				$median_spacer_length=$arr_spacer_lengths[int($middle_pos-1)];
			}		
				
		return ($median_spacer_length);
				
	}	
	

sub get_all_cds_or_crispr_positions()
	{
		my($accession,$all_gene_positions_file,$search_term,$return_array)=@_;
		
		if($all_gene_positions_file eq "NA"){return 1;}
		
		#---- first check if the accession is present in the file or not
		my @arr_gene_positions;
		
		#print "\ngrep -w '$accession' $tmp_dir\/$all_gene_positions_file | grep -w '$search_term'\n\n";
		@arr_gene_positions=`grep -w '$accession' $tmp_dir\/$all_gene_positions_file | grep -w '$search_term' >&1`;
		
		if($#arr_gene_positions<=0)
			{
				#@arr_gene_positions=`grep -w '$search_term' $tmp_dir\/$all_gene_positions_file >&1`;
			}	
		
		foreach my $line(@arr_gene_positions)
			{
				push(@{$return_array},$line);
			}
			
		return 1;	
	}

sub check_array_existance()
	{
				my($array_quality_score,$range_start,$range_stop,$tmp_output_file)=@_;
				
				($range_start,$range_stop)=&sort_two_numbers($range_start,$range_stop);
				my $current_middle_point=int($range_start+($range_stop-$range_start)/2);
				
				my $already_exists=0;
				my $existing_q_score=0;
				
				open(RD5,"$tmp_dir\/$tmp_output_file") or print "$!";
				flock(RD5,2);
				my @arr_rd5=<RD5>;
				close(RD5);
				
				#foreach my $line(@arr_rd5)
				my @arr_lines_to_delete;
				my $lines_to_be_deleted=0;
				for(my $i=0;$i<=$#arr_rd5;$i++)
					{
						$line=$arr_rd5[$i];
						chomp $line;$line=~s/\r//g; if(not $line or $line eq ""){next;}
						my @arr_line=split('%%%',$line);
						my $accession=$arr_line[0];
						my $start=$arr_line[1];
						my $stop=$arr_line[2];
						my $q_score=$arr_line[3];
						my $existing_array=$arr_line[4];
						
						($start,$stop)=&sort_two_numbers($start,$stop);
						my $existing_middle_point=int($start+($stop-$start)/2);
						
						#if(($current_middle_point>$start and $current_middle_point<$stop) or ($existing_middle_point>$range_start and $existing_middle_point<$range_stop))
						if(($range_start>=$start and $range_start<=$stop) or ($range_stop>=$start and $range_stop<=$stop) or ($start>=$range_start and $start <= $range_stop) or ($stop>=$range_start and $stop<=$range_stop) )
							{
								#--- an array exist for this range
								if($array_quality_score<=$q_score)
									{
										$already_exists=1;
										$existing_q_score=$q_score;
										last;
									}
								else{
										#--- the current version is better than existing one, so delete the row------
										#my $j=$i+1;
										#push(@arr_lines_to_delete,$j);
										$lines_to_be_deleted=1;
										$arr_rd5[$i]="";
									}
							}
						
					}
					
				#---- now write the remaining 	
				if($lines_to_be_deleted!=0)
					{
						open(WR5,">$tmp_dir\/$tmp_output_file") or print "$!";
						flock(WR5,2);
						for(my $i=0;$i<=$#arr_rd5;$i++)
							{
								if(defined $arr_rd5[$i] and $arr_rd5[$i] ne "")
									{
										$line=$arr_rd5[$i];
										print WR5 "$line\n";
									}	
							}
						close(WR5);		
								
					}	
					
				return($already_exists,$existing_q_score);
	}
	

sub run_water()
	{
		my($range,$accession,$seq1,$seq2,$gapopen,$gapextend)=@_;
		
		
		my $time_1 = &get_unique_id();
				
		my $outfile=$time_1."_".$accession.$range."_output.txt";
		open(WR,">$tmp_dir\/$outfile");close(WR);	#if(-e "$tmp_dir\/$outfile"){system("chmod 777 $tmp_dir\/$outfile");}
		
		if(not defined $seq1 or $seq1 eq "" or $seq1=~/[^ACGTU]/i){return($outfile);}
		if(not defined $seq2 or $seq2 eq "" or $seq2=~/[^ACGTU]/i){return($outfile);}

		

						
		my $file_seq_1=$time_1."_".$accession.$range."_spacer_bases.txt";
		my $file_seq_2=$time_1."_".$accession.$range."_repeat_bases.txt";
		
		
		open(WR,">$tmp_dir\/$file_seq_1");flock(WR,2);print WR "$seq1";close(WR); #if(-e "$tmp_dir\/$file_seq_1"){system("chmod 777 $tmp_dir\/$file_seq_1");}
		open(WR,">$tmp_dir\/$file_seq_2");flock(WR,2);print WR "$seq2";close(WR); #if(-e "$tmp_dir\/$file_seq_2"){system("chmod 777 $tmp_dir\/$file_seq_2");}
		
		
		if(-e "$tmp_dir\/$file_seq_1" and -e "$tmp_dir\/$file_seq_2")
			{				
				#system("chmod 777 $tmp_dir\/$file_seq_1");
				#system("chmod 777 $tmp_dir\/$file_seq_2");
				
				my $ret_msg=`water -asequence $tmp_dir\/$file_seq_1 -bsequence $tmp_dir\/$file_seq_2 -gapopen $gapopen -gapextend $gapextend -outfile $tmp_dir\/$outfile >&1 2>&1`;
						#print "./water -asequence $tmp_dir\/$file_seq_1 -bsequence $tmp_dir\/$file_seq_2 -gapopen $gapopen -gapextend $gapextend -outfile $tmp_dir\/$outfile >/dev/null 2>&1<br>";
				if($ret_msg=~/core dumped/i)
					{
						#print "\$seq1=$seq1\n\$seq2=$seq2\n";
					}	
					
				#print "\n\$seq1=$seq1\n\$seq2=$seq2\n\n";		
			}
		else{
				sleep(1);
				
				open(WR,">$tmp_dir\/$outfile");close(WR);
				open(WR,">$tmp_dir\/$file_seq_1");flock(WR,2);print WR "$seq1";close(WR); #if(-e "$tmp_dir\/$file_seq_1"){system("chmod 777 $tmp_dir\/$file_seq_1");}
				open(WR,">$tmp_dir\/$file_seq_2");flock(WR,2);print WR "$seq2";close(WR); #if(-e "$tmp_dir\/$file_seq_2"){system("chmod 777 $tmp_dir\/$file_seq_2");}
				
				my $ret_msg=`water -asequence $tmp_dir\/$file_seq_1 -bsequence $tmp_dir\/$file_seq_2 -gapopen $gapopen -gapextend $gapextend -outfile $tmp_dir\/$outfile >&1 2>&1`;
						#print "./water -asequence $tmp_dir\/$file_seq_1 -bsequence $tmp_dir\/$file_seq_2 -gapopen $gapopen -gapextend $gapextend -outfile $tmp_dir\/$outfile >/dev/null 2>&1<br>";
				if($ret_msg=~/core dumped/i)
					{
						print "\$seq1=$seq1\n\$seq2=$seq2\n";
					}
			}	
			

							
		#system("echo '$ret_msg' >>log1.txt");
						
		unlink("$tmp_dir\/$file_seq_1");					
		unlink("$tmp_dir\/$file_seq_2");
		
		return($outfile);
	}
sub print_array()
	{
		my($model_repeat,$current_array)=@_;
		foreach my $l(@{$current_array})
			{
				print "$l\n";
			} 
		print "\t$model_repeat\n\n";	
	}	
	
#-- end ---------------------------------
return 1;


