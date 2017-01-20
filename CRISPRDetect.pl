#!/usr/bin/perl
our $version=2.2;

#------------------------------------------------------------------------------------------------------------------------------------------------
#Program : 	CRISPRDetect_v2.2
#Author: 	Ambarish Biswas
#Contact:	Chris Brown [chris.brown@otago.ac.nz] or Ambarish Biswas [ambarish.biswas@otago.ac.nz]
#------------------------------------------------------------------------------------------------------------------------------------------------
###### CRISPRDetect_v2.1 supports a range of parameters passed through command line without changing this script. 
###### The parameters are named in relation to the function/module it is related with. For example, the parameter 
###### 		'extend_array' is associated with extending the CRISPR arrays, which uses 'ea_allowed_percent_similarity' and 
###### 		'ea_dynamic_search'. Note the 'ea_' in the beginning which is abbreviated from 'extend array'. The parameters are
###### 		listed in the order of the modules listed in the CRISPRDetect webserver, and can be easily understood by comparing with the website.
#------------------------------------------------------------------------------------------------------------------------------------------------
use Cwd 'abs_path';
our $cd_path=abs_path($0); $cd_path=&get_path($cd_path);
require "$cd_path/CD_MODULES/CRISPRDETECT_SUBS_1.pm";
use lib '$cd_path/lib';


use strict;
use Term::ANSIColor;





our $tmp_dir="$cd_path/tmp";
my $no_of_threads=4;


use Parallel::ForkManager;
my $pm = new Parallel::ForkManager($no_of_threads);      # specify the number of parallel processes CRISPRDetect should use. Defaults to 4.










#----- before running check all dependencies are installed -----------------------------------------------------------------------------

my($clustalw_found,$water_found,$seqret_found,$RNAfold_found,$cd_hit_found,$blastn_found)=&check_executables_required_for_CRISPRDetect(); #----- check the required 3rd party executables ------
if($clustalw_found=~/^Not found/ or $water_found=~/^Not found/ or $seqret_found=~/^Not found/ or $RNAfold_found=~/^Not found/ or $cd_hit_found=~/^Not found/ or $blastn_found=~/^Not found/)
	{
		print "\nError: Some of the dependencies not found in system path. Program terminating.\n"; exit;
	}
if(not -e $tmp_dir)
	{		
		system("mkdir -p $tmp_dir");system("chmod 775 $tmp_dir");
	}












#******************************************************** change the parameters for initial (putative) CRISPR prediction *****************************************
my $keep_input_sequences=0;

my $word_length=11;
my $minimum_word_repeatation=3;

my $minimum_spacer_gap=50-$word_length;
my $maximum_spacer_gap=100+$word_length;


my $max_gap_between_crisprs=500;

my $left_flank_length=500;
my $right_flank_length=500;

my $continue_from_last_process=0;
my $continue_from_last_process_file="already_processed_sequences.list";



#-------------------------------- filtering option
my $repeat_length_cutoff=17;
my $minimum_no_of_repeats=3;
my $minimum_repeat_length=23;
my $array_quality_score_cutoff=4;

#-------------------------------- first turn on/off specific modules auto execution using these control parameters

my $remove_insertion_from_repeat=1;

my $fix_gaps_in_repeats=1;

my $extend_array=1;
	my $ea_allowed_percent_similarity=67;
	my $ea_dynamic_search=1;

my $search_unidentified_repeat_in_spacer_sequence=1;
	my $allowed_percent_similarity=67;	
	my $su_dynamic_search=0;
	
my $trim_repeats=1;
	my $trimming_cutoff="AUTO";
	my $user_side_to_trim="NA";
	my $user_no_of_bases_to_trim=0;


my $increase_repeat_length=1;
	my $repeat_extension_identity="AUTO";    		#-- default >50% --------------------
	my $user_side_to_increase_length="NA";
	my $user_no_of_bases_to_increase_length=0;
	
	
my $shorten_array=1;
	my $user_side_to_shorten="NA";
	my $sa_allowed_percent_similarity=66;


my $check_consensus=1;						# ---by default set to 1


my $find_alternate_repeat=1;				# ---default set to 1 
	my $potential_alternate_repeat="NA";


my $check_direction=1;						# ---default set to 1 

	#-------------------------- change the parameters for directional analysis directional analysis (CRISPRDirection) ----------------------------------------------------------------------------

	my $check_motif_in_repeat=1; 
		my $Motif_match_score=4.50;			# Default: Score 4.5 (acts as a fileter); If you choose to use this method like the other methods make the score 0.50 [i.e. PPV (1) - 0.50];
		my $motif="ATTGAAA.?";				# Default: "ATTGAAA.?" where '.?' is equal to N, another example is "A.?.?TGAAA[C|G]" which is same as ANNTGAAAC or ANNTGAAAG 

	my $check_A_and_T_ratio_in_repeat=1;
		my $A_and_T_ratio_score=0.37;		# Default: PPV (0.87) - 0.50
		
	my $check_similarity_with_reference_repeat=1; 
		my $Similarity_score=4.50;			# Default: Score 0.50 [i.e. PPV (1) - 0.50]; If you want to make it act as a filter use Score 4.5 or higher; 	
		my $allowed_no_of_mismatches=3;		# Default 3		#Note: We noticed upto allowed 6 bases mismatch there were no false predictions (Refer paper/supplement)
		
	my $check_secondary_structure_of_repeat=1;		
		my $MFE_score=0.37;					# Default: PPV (0.87) - 0.50		
		my $MFE_cutoff=1;					# Default 1
		my $MFE_minimum_difference=0.25;	# Default 0.25
		my $MFE_exclude_bases=5;			# Default 5	

	my $check_array_degeneracy=1;	 
		my $array_degeneracy_score=0.41;	# Default: PPV (0.91) - 0.50
		my $permitted_mutation_per_array=0;	# Default 0

	my $check_AT_distribution_in_flanks=1; 
		my $AT_distribution_score=0.27;		# Default: PPV (0.77) - 0.50
		my $AT_distribution_window=60;		# Default: 60
		my $AT_distribution_minimum_percentage_difference=10;	# Default: 10 


	my $check_longer_leader=1;
		my $Longer_leader_score=0.18;		# Default: PPV (0.68) - 0.50



my $user_reverse_the_array=0;					# ---default set to 0

###################################################################### end of optional user input section #####################################################


































##################################################################### do not change anyting below this, unless you know what you doing #######################################################







#--------- open library of repeats with confirmed direction -----------------------
my %lib_of_repeats_with_confirmed_direction;				
open(LIB_REP,"$cd_path/Ref_lib_files/verified_repeats_with_family.txt") or print "$!";
my @arr_ref_lib=<LIB_REP>;
close(LIB_REP);
				
foreach my $ref_line(@arr_ref_lib)
	{
		if($ref_line=~/^#/){next;}  # skip the comment lines if any
		
		chomp $ref_line;$ref_line=~s/\r//g;
		my($rep,$fam_type)=split('\t',$ref_line);
		my $translated_rep=$rep; $translated_rep=~tr/T/U/;
		$lib_of_repeats_with_confirmed_direction{$rep}=$fam_type;
		$lib_of_repeats_with_confirmed_direction{$translated_rep}=$fam_type;
	}


my $blast_db_file_of_known_repeats="$cd_path/BLAST_DB_OF_KNOWN_REPEATS/known_repeats.fa";





#---------------------------------------------------------------------------------------

#my $input_array_file;

my $all_gene_positions_folder="$tmp_dir\/";
#my $all_gene_positions_file="NA";

my $output_file="$tmp_dir\/CRISPRDetect_output_".&get_unique_id().".txt";
my $ref_lib_file="NA";
my $gff_file="NA";
my $filtered_out_crisprs="NA";















my $quiet=0;

my $species="Bacteria";


my $idnumber=&get_unique_id();

my $all_gene_positions_file="NA";#$idnumber."_CDS_positions.txt"; open(WR,">$tmp_dir\/$all_gene_positions_file");close(WR);

my $input_sequence_file=$idnumber."_fasta_sequence.txt";
#my $output_file;
my @arr_sequence_files;
my %hash_id_lookup_table;
my %hash_of_accession_and_species;


my $repeats_file=$idnumber."_repeats.txt";;
my $process_repeats_only=0;

my $array_seq_range="NA";

#------ step 0: get the input file ----------------------------------------
if(not defined $ARGV[0]){&show_help($clustalw_found,$water_found,$seqret_found,$RNAfold_found,$cd_hit_found,$blastn_found);exit;}

my $output_file_specified=0;

for(my $i=0;$i<=$#ARGV;$i++)
	{
		#------ silencing the program's step by step reporting
		if($ARGV[$i]=~/-tmp_dir/)
			{
				$tmp_dir=$ARGV[$i+1];
				if(not -e $tmp_dir)
					{
						print "\nError: The specified directory: $tmp_dir not found. User should create the directory, and give read and write permission to perl and its dependent tools. Best to use the default '/tmp' directory.\n\n";exit;
					}	
												
			}
		elsif($ARGV[$i]=~/-keep_input_sequences/)
			{
				$keep_input_sequences=$ARGV[$i+1];	
				
				$keep_input_sequences=abs(int($keep_input_sequences));
								
			}
		elsif($ARGV[$i]=~/-quiet$/)
			{
				$quiet=$ARGV[$i+1];	
				
				$quiet=abs(int($quiet));
								
			}
		elsif($ARGV[$i]=~/-q$/)
			{
				$quiet=$ARGV[$i+1];	
				
				$quiet=abs(int($quiet));
								
			}
		#elsif($ARGV[$i]=~/-all_gene_positions_folder/)
		#		{
		#			$all_gene_positions_folder=$ARGV[$i+1];					
		#		}
		elsif($ARGV[$i]=~/-all_gene_positions_file/)
				{
					$all_gene_positions_file=$ARGV[$i+1];					
				}
		elsif($ARGV[$i]=~/-ref_lib_file/)
				{
					$ref_lib_file=$ARGV[$i+1];	
					chomp $ref_lib_file; $ref_lib_file=~s/\r//g;	
					
					#print "Found ref_lib_file\n\n";	
					
					if(-e $ref_lib_file and $ref_lib_file!~/NA/)
						{
							#print "\n\n\$ref_lib_file= $ref_lib_file\n\n";
							open(LIB_REP,"$ref_lib_file") or print "$!";
							my @arr_ref_lib=<LIB_REP>;
							close(LIB_REP);
							
							#--- now create a blast_db_of_known_repeats.fa file	
							my $tmp_blast_db_file_of_known_repeats="$tmp_dir\/".&get_unique_id()."_known_repeats.fa";	
							
							system("cat $blast_db_file_of_known_repeats >$tmp_blast_db_file_of_known_repeats");
							
							open(UR,">>$tmp_blast_db_file_of_known_repeats");					
							my $user_repeat_index=0;
							foreach my $ref_line(@arr_ref_lib)
								{
									if($ref_line=~/^#/){next;}  # skip the comment lines if any
									$user_repeat_index++;
									
									
									chomp $ref_line;$ref_line=~s/\r//g;
									
									my($rep,$fam_type)=split('\t',$ref_line);
									#print "$ref_line: [$rep,$fam_type]\n";
									my $translated_rep=$rep; $translated_rep=~tr/T/U/;
									$lib_of_repeats_with_confirmed_direction{$rep}=$fam_type;
									$lib_of_repeats_with_confirmed_direction{$translated_rep}=$fam_type;
									
									if($fam_type=~/\S{1,}/)
										{
											$fam_type=~s/_/-/g;
											print UR ">UR_$user_repeat_index\_$fam_type\n$rep\n";
											#print ">UR_$user_repeat_index\_$fam_type\n$rep\n";
										}
									else{	
											print UR ">UR_$user_repeat_index\_USER-REPEAT\n$rep\n";
										}	
								}
							close(UR);			
							
							
							#--- now make a blast DB
							system("makeblastdb -in $tmp_blast_db_file_of_known_repeats -dbtype nucl -hash_index >/dev/null 2>/dev/null");
								
							#---- now update the $blast_db_file_of_known_repeats
							$blast_db_file_of_known_repeats=$tmp_blast_db_file_of_known_repeats;	
							
							#print "\nNew \$blast_db_file_of_known_repeats= $blast_db_file_of_known_repeats\n\n";
						}
					
						
				}		
		#------ set the no. of parallel processes
		elsif($ARGV[$i]=~/-T$/)
			{
				$no_of_threads=$ARGV[$i+1];
				$no_of_threads=int($no_of_threads);	
				
				if(-e "/proc/cpuinfo")
					{				
						my $no_of_system_cpus=`grep 'processor' /proc/cpuinfo | wc -l >&1`;
						chomp $no_of_system_cpus; $no_of_system_cpus=~s/\r//g;
						$no_of_system_cpus=int($no_of_system_cpus);
						
						if($no_of_system_cpus >1)
						{
							if($no_of_threads < $no_of_system_cpus)
								{
									#-- sweet
									if($no_of_threads ==0)
										{
											$no_of_threads=$no_of_system_cpus-1;
										}
								}
							elsif($no_of_threads>=$no_of_system_cpus)
								{
									$no_of_threads=$no_of_system_cpus-1;								
								}
							else{
									$no_of_threads=1;
								}	
						}	
					}
				else{
						$no_of_threads=1;
					}	
				
				#--- now reinitiate the $pm
				$pm = new Parallel::ForkManager($no_of_threads);		
								
			}
				
			
		#--------------------------- search CRISPR-Hotspots ---------------------------------------------------
		elsif($ARGV[$i]=~/-word_length$/)
			{
				$word_length=$ARGV[$i+1];	
				
				$word_length=abs(int($word_length));
				
				if($word_length<6){$word_length=6;}
				
				$minimum_spacer_gap=50-$word_length;
				$maximum_spacer_gap=100+$word_length;				
			}
		elsif($ARGV[$i]=~/-minimum_word_repeatation$/)
			{
				$minimum_word_repeatation=$ARGV[$i+1];					
			}		
		#--------------------------- input options ------------------------------------------------------------
		elsif($ARGV[$i]=~/-minimum_no_of_repeats/)
			{
				$minimum_no_of_repeats=$ARGV[$i+1];					
			}
		elsif($ARGV[$i]=~/-repeat_length_cutoff/)
			{
				$repeat_length_cutoff=$ARGV[$i+1];					
			}	
			
		elsif($ARGV[$i]=~/-g$/)
			{
				#----- input is a gbk file : first extract sequence from it------------------------------------
				my $filename=$ARGV[$i+1];			
				
				if($filename=~/\.gbk$/)
					{
						my($accession,$defination,$species)=&extract_information_and_sequence_from_gbk_file($filename);							

						push(@arr_sequence_files,$accession);
						$hash_id_lookup_table{$accession}="$accession-$defination";	
						$hash_of_accession_and_species{$accession}=$species;
						##------------ get the accession ---------------------------------------------------------------
						#&process_fasta_file($tmp_file,\@arr_sequence_files,\%hash_id_lookup_table);		
						
								
						#----------- now get the gene positions -------------------------------------------------------
						if($all_gene_positions_file eq "NA")
							{
								$all_gene_positions_file=&get_unique_id()."_CDS_positions.txt";							
							}			
						&process_gbk_file($filename,$all_gene_positions_file);
						#---------------------------------------------------------------------------------------------- 
					}
				elsif($filename=~/\.gb$/)
					{
						my($accession,$defination,$species)=&extract_information_and_sequence_from_gbk_file($filename);							

						push(@arr_sequence_files,$accession);
						$hash_id_lookup_table{$accession}="$accession-$defination";	
						$hash_of_accession_and_species{$accession}=$species;
						##------------ get the accession ---------------------------------------------------------------
						#&process_fasta_file($tmp_file,\@arr_sequence_files,\%hash_id_lookup_table);		
						
								
						#----------- now get the gene positions -------------------------------------------------------
						if($all_gene_positions_file eq "NA")
							{
								$all_gene_positions_file=&get_unique_id()."_CDS_positions.txt";							
							}			
						&process_gbk_file($filename,$all_gene_positions_file);
						#----------------------------------------------------------------------------------------------
					}	
				elsif($filename=~/\.gbff$/)
					{
						#----split the gbff file to multiple contigs file
						my @arr_gbk_files;
						&split_gbff_file_to_individual_gbk_files($filename,\@arr_gbk_files);
						
						foreach my $gbk_file(@arr_gbk_files)
							{
								#print "$gbk_file\n";
								my($accession,$defination,$species)=&extract_information_and_sequence_from_gbk_file($gbk_file);								

								push(@arr_sequence_files,$accession);
								
								$hash_id_lookup_table{$accession}="$accession-$defination";	
								$hash_of_accession_and_species{$accession}=$species;
								##------------ get the accession ---------------------------------------------------------------
								#&process_fasta_file($tmp_file,\@arr_sequence_files,\%hash_id_lookup_table);		
								
										
								#----------- now get the gene positions -------------------------------------------------------
								if($all_gene_positions_file eq "NA")
									{
										$all_gene_positions_file=&get_unique_id()."_CDS_positions.txt";							
									}			
								&process_gbk_file($gbk_file,$all_gene_positions_file);
								#----------------------------------------------------------------------------------------------
								
								unlink($gbk_file);
							}
							
					}		

			}
		elsif($ARGV[$i]=~/-f$/)
			{
				#----- input is a fasta file ------------------------------------------------------------------
				my $tmp_file=$ARGV[$i+1];
				
				#--- now copy the sequence from $source_fasta_file to $tmp_dir\/$input_sequence_file
				&process_fasta_file($tmp_file,\@arr_sequence_files,\%hash_id_lookup_table);
			}
		elsif($ARGV[$i]=~/-r$/)
			{
				#----- input is a directory of GBK files  ---------------------------------------------------
				my $gbk_folder=$ARGV[$i+1];
				
				my @all_gbk_files=`find $gbk_folder -name '*.gbk' >&1`;
				
				my $total_gbk_files_to_process=$#all_gbk_files+1;
				print "Total $#all_gbk_files to process...\n";
				
				#---- write $accession,$defination,$species in a tmp file, and load up later (as child prosses cant write to parent array or hashes)
				my $tmp_acc_def_and_species_file=&get_unique_id()."_acc_def_species_list.txt";
				system("echo '' >$tmp_dir\/$tmp_acc_def_and_species_file");
				
				if($all_gene_positions_file eq "NA")
					{
						$all_gene_positions_file=&get_unique_id()."_CDS_positions.txt";							
					}
				
				my $file_count=0;
				foreach my $filename(@all_gbk_files)
					{
						$file_count++;
						#if($file_count>25){last;}
						#if($filename!~/NC_000909/){next;}
						
						#if($file_count % 100 == 0)
						#	{
						#		sleep(int(rand(10)));
						#	}
						
						chomp $filename;$filename=~s/\r//g;
						
						$total_gbk_files_to_process--;
						select(undef, undef, undef, 0.25); #--- will sleep for 1/4 seconds
						$pm->start and next; # do the fork
						print "\rRemaining: $total_gbk_files_to_process\t processing $filename...\t\t\t\r";
						
						my($accession,$defination,$species);
						$accession="NA";
						
						my $no_of_tries=0;
						while($accession =~ /NA/)
							{
								$no_of_tries++;
								
								($accession,$defination,$species)=&extract_information_and_sequence_from_gbk_file($filename);				
								
								if($no_of_tries>10){last;}
							}	
				
						if($accession ne "NA")
							{
								#push(@arr_sequence_files,$accession);
								#$hash_id_lookup_table{$accession}="$accession-$defination";	
								#$hash_of_accession_and_species{$accession}=$species;
							
								open(APP,">>$tmp_dir\/$tmp_acc_def_and_species_file") or print "$!";	
								flock(APP,2);
								print APP "$accession\t$defination\t$species\n";
								close(APP);	
								#----------- now get the gene positions -------------------------------------------------------
											
								&process_gbk_file($filename,$all_gene_positions_file);
								#----------------------------------------------------------------------------------------------
						
								
							}						
						
						#print "$filename\n";
						
						$pm->finish;
					}
				 $pm->wait_all_children;
				 
				 #print "\n";
				 
				 
				 #------- now load up the acc and other details ---
				 open(RD,"$tmp_dir\/$tmp_acc_def_and_species_file") or print "$!";	
				 flock(RD,2);
				 while(my $line=<RD>)
					{
						chomp $line;$line=~s/\r//g; if(not $line){next;}
						my($accession,$defination,$species)=split('\t',$line);
						
						push(@arr_sequence_files,$accession);
						$hash_id_lookup_table{$accession}="$accession-$defination";	
						$hash_of_accession_and_species{$accession}=$species;
						
					}				 
				 close(RD);
			}	
		
		elsif($ARGV[$i]=~/-continue$/)
			{
				#print "Continue set to true\n\n";
				#----- input is a fasta file ------------------------------------------------------------------
				$continue_from_last_process=$ARGV[$i+1];
				if($continue_from_last_process==0 or not -e "$tmp_dir\/$continue_from_last_process_file")
					{
						open(CLEAR,">$tmp_dir\/$continue_from_last_process_file") or print "$!\n";close(CLEAR);
					}
			}
		#--------- output options -----------------------------------------------------------------------------
		elsif($ARGV[$i]=~/-o$/)
			{
				#----- input is a fasta repeat(s) only file ---------------------------------------------------
				$output_file=$ARGV[$i+1];
				
				
				if(not defined $ARGV[$i+1])
					{
						print "\nWrong input. A filename is expected after -o [e.g. -o output.txt]\n\n";
						exit;
					}
				else{
						open(WR,">$output_file") or print "$!";close(WR);system("chmod 777 $output_file");
						
						$gff_file=$output_file.".gff";
						open(WR,">$gff_file") or print "$!";close(WR);system("chmod 777 $gff_file");
						
						$filtered_out_crisprs=$output_file.".fp";
						open(WR,">$filtered_out_crisprs") or print "$!";close(WR);system("chmod 777 $filtered_out_crisprs");
					}	
				$output_file_specified=1;	
			}
		#----------------------------------------------------------			
	
			#-----------------------------------------------------------	
			elsif($ARGV[$i]=~/-species/)
				{
					$species=$ARGV[$i+1];					
				}	
				
			elsif($ARGV[$i]=~/-left_flank_length/)
				{
					$left_flank_length=$ARGV[$i+1];					
				}
			elsif($ARGV[$i]=~/-right_flank_length/)
				{
					$right_flank_length=$ARGV[$i+1];					
				}	
			elsif($ARGV[$i]=~/-max_gap_between_crisprs/)
				{
					$max_gap_between_crisprs=$ARGV[$i+1];					
				}				
				
			#--------------	
			elsif($ARGV[$i]=~/-remove_insertion_from_repeat/)
				{
					$remove_insertion_from_repeat=$ARGV[$i+1];					
				}
				
			elsif($ARGV[$i]=~/-extend_array/)
				{
					$extend_array=$ARGV[$i+1];					
				}
			elsif($ARGV[$i]=~/-ea_allowed_percent_similarity/)
				{
					$ea_allowed_percent_similarity=$ARGV[$i+1];					
				}
			elsif($ARGV[$i]=~/-ea_dynamic_search/)
				{
					$ea_dynamic_search=$ARGV[$i+1];					
				}
					
			elsif($ARGV[$i]=~/-fix_gaps_in_repeats/)
				{
					$fix_gaps_in_repeats=$ARGV[$i+1];					
				}
				
			elsif($ARGV[$i]=~/-search_unidentified_repeat_in_spacer_sequence/)
				{
					$search_unidentified_repeat_in_spacer_sequence=$ARGV[$i+1];					
				}
			elsif($ARGV[$i]=~/-su_dynamic_search/)
				{
					$su_dynamic_search=$ARGV[$i+1];					
				}	
				
			elsif($ARGV[$i]=~/-allowed_percent_similarity/)
				{
					$allowed_percent_similarity=$ARGV[$i+1];					
				}				
			elsif($ARGV[$i]=~/-trim_repeats/)
				{
					$trim_repeats=$ARGV[$i+1];					
				}
			elsif($ARGV[$i]=~/-minimum_repeat_length/)
				{
					$minimum_repeat_length=$ARGV[$i+1];					
				}	
				
			elsif($ARGV[$i]=~/-trimming_cutoff/)
				{
					$trimming_cutoff=$ARGV[$i+1];					
				}
				
			elsif($ARGV[$i]=~/-user_side_to_trim/)
				{
					$user_side_to_trim=$ARGV[$i+1];					
				}
				
			#elsif($ARGV[$i]=~/-user_no_of_bases_to_trim/)
			#	{
			#		$user_no_of_bases_to_trim=$ARGV[$i+1];					
			#	}
				
			elsif($ARGV[$i]=~/-increase_repeat_length/)
				{
					$increase_repeat_length=$ARGV[$i+1];					
				}
			
			elsif($ARGV[$i]=~/-repeat_extension_identity/)
				{
					$repeat_extension_identity=$ARGV[$i+1];					
				}	
				
			elsif($ARGV[$i]=~/-user_side_to_increase_length/)
				{
					$user_side_to_increase_length=$ARGV[$i+1];					
				}
				
			elsif($ARGV[$i]=~/-user_no_of_bases_to_increase_length/)
				{
					$user_no_of_bases_to_increase_length=$ARGV[$i+1];					
				}
				
				
			#elsif($ARGV[$i]=~/-check_consensus/)
			#	{
			#		$check_consensus=$ARGV[$i+1];					
			#	}
							
			#elsif($ARGV[$i]=~/-find_alternate_repeat/)
			#	{
			#		$find_alternate_repeat=$ARGV[$i+1];					
			#	}
							
			elsif($ARGV[$i]=~/-check_direction/)
				{
					$check_direction=$ARGV[$i+1];					
				}
			#elsif($ARGV[$i]=~/-ref_lib_file/)
			#	{
			#		$ref_lib_file=$ARGV[$i+1];	
			#		$ref_lib_file="$tmp_dir\/$ref_lib_file";				
			#	}	




			elsif($ARGV[$i]=~/-shorten_array/)
				{
					$shorten_array=$ARGV[$i+1];					
				}	
				
			elsif($ARGV[$i]=~/-sa_allowed_percent_similarity/)
				{
					$sa_allowed_percent_similarity=$ARGV[$i+1];					
				}				
					
			elsif($ARGV[$i]=~/-user_side_to_shorten/)
				{
					$user_side_to_shorten=$ARGV[$i+1];					
				}
			elsif($ARGV[$i]=~/-array_quality_score_cutoff/)
				{
					$array_quality_score_cutoff=$ARGV[$i+1];					
				}	
		#----- optional, will be removed later ----------------------------------------------------------------
		#elsif($ARGV[$i]=~/-array_seq_range/)
		#	{
		#		#----- input is a fasta repeat(s) only file ---------------------------------------------------
		#		$array_seq_range=$ARGV[$i+1];
		#		
		#	}	
		elsif($ARGV[$i]=~/-h/ or $ARGV[$i]=~/-help/)
			{
				#----- CRISPRDetect help ---------------------------------------------------
				&show_help($clustalw_found,$water_found,$seqret_found,$RNAfold_found,$cd_hit_found,$blastn_found);
				
				exit;
			}	
	}


if($output_file_specified==0)
	{
		print "\nWrong input. A output filename should be specified [e.g. -o output.txt]\n\n";
		exit;
	}



###################################################################




#----- step 1: predict CRISPR hotspots ---------------------



my $date_n_time=localtime(time);
my $file_index=0;
my $combined_hotspots_file=$idnumber."_combined_hotspots.txt";
system("echo '# Thank you for using CRISPRDetect' >$tmp_dir\/$combined_hotspots_file"); 
system("echo '# Author: Ambarish Biswas' >$tmp_dir\/$combined_hotspots_file"); 
system("echo '# $date_n_time' >$tmp_dir\/$combined_hotspots_file");
system("echo '' >$tmp_dir\/$combined_hotspots_file");
system("echo '' >$tmp_dir\/$combined_hotspots_file");
system("chmod 777 $tmp_dir\/$combined_hotspots_file");



my $remaining_sequences=$#arr_sequence_files+1;

#print "$combined_hotspots_file\n\n";
foreach my $seq_file(@arr_sequence_files)
	{		
			$remaining_sequences--;
			$file_index++;
			#print "\n$remaining_sequences\t sequences remaining for finding putative arrays\n ";
			
			if(-e "$tmp_dir\/$combined_hotspots_file")
				{
					my $hotspots_identified=`grep '>' $tmp_dir\/$combined_hotspots_file | wc -l >&1`;
					$hotspots_identified=int($hotspots_identified);
					
					if($hotspots_identified>0)
						{
							if($quiet !=1)
								{
									print "\r$hotspots_identified putative CRISPRs identified";
								}	
						}					
				}
			
			#------ check if the sequence already been processed b4, if continue is set
			if(-e "$tmp_dir\/$continue_from_last_process_file" and $continue_from_last_process==1)
				{
					#print "grep -w '$seq_file' $tmp_dir\/$continue_from_last_process_file\n";
					my $is_already_processed=`grep -w '$seq_file' $tmp_dir\/$continue_from_last_process_file >&1`;
					if($is_already_processed=~/$seq_file/)
						{
							print "$seq_file already processed. Continuing..\n";
							next;
						}	
				}	
			
			#print "$seq_file\n";
			select(undef, undef, undef, 0.25); #--- will sleep for 1/4 seconds
			$pm->start and next;
			
			$seq_file=$seq_file."\.fna";
			
			#&predict_crispr_hotspots($seq_file,$combined_hotspots_file);	
			
			my ($total_hotspots,$potential_crisprs)=&get_crisprdetect_hotspots($word_length,$minimum_word_repeatation,$minimum_spacer_gap,$maximum_spacer_gap,$seq_file,$combined_hotspots_file,\%hash_id_lookup_table);	
			
			#print "Total $potential_crisprs potential CRISPR hotspots identified out of $total_hotspots hotspots\n";
			
			$pm->finish;
	}
	
	$pm->wait_all_children;

print "\n";
#exit;


#foreach my $acc(keys %hash_id_lookup_table)
#	{
#		#print "$acc\n";
#	}







#----- step 3: Process the CRISPR arrays to get the longest one and load up %hash_of_all_crispr_hotspots -----------
my %hash_of_all_crispr_hotspots;

#print "\nGoing to load up the CRISPR hotspots ...";

my $total_records=&load_crispr_hotspots("$tmp_dir\/$combined_hotspots_file",\%hash_of_all_crispr_hotspots);

unlink("$tmp_dir\/$combined_hotspots_file");
#print ".. Done.\n\n";	




#exit;		
		
	

				


#------- get the total CRISPR hotspots
foreach my $accession(sort keys %hash_of_all_crispr_hotspots)
	{
		foreach my $range(sort keys %{$hash_of_all_crispr_hotspots{$accession}})
			{
				#print "$accession-$range\n";
				$total_records++;
			}	
	}
#exit;
#----- step 4: pass the hash of CRISPR arrays to check_orientation ---------


#--- get all the processed accessions from output file ----------------------

#my @already_processed_arrays=`grep '>' compiled_CRISPRDetect_output.txt | awk '{FS="-";print \$1}' | awk '{FS=">";print \$2}' |  awk '{arr[\$NF] = \$0} END { for (key in arr) { print arr[key] } }' >&1`;
#my %hash_of_already_processed_arrays;
#foreach my $accession(@already_processed_arrays)
#	{
#		chomp $accession;$accession=~s/\r//g;
#		$hash_of_already_processed_arrays{$accession}=1;
#	}
#my @already_processed_arrays2=`grep '>' compiled_CRISPRDetect_output.txt.fp | awk '{FS="-";print \$1}' | awk '{FS=">";print \$2}' |  awk '{arr[\$NF] = \$0} END { for (key in arr) { print arr[key] } }' >&1`;
#foreach my $accession(@already_processed_arrays2)
#	{
#		chomp $accession;$accession=~s/\r//g;
#		$hash_of_already_processed_arrays{$accession}=1;
#	}

if($quiet !=1)
	{
		print "\tTotal CRISPRs to process: $total_records\n";
	}
foreach my $accession(sort keys %hash_of_all_crispr_hotspots)
 {	
	
	
	#------ check if the sequence already been processed b4, if continue is set ------------------------
	if(-e "$tmp_dir\/$continue_from_last_process_file" and $continue_from_last_process==1)
		{
			#print "grep -w '$accession' $tmp_dir\/$continue_from_last_process_file\n";
			my $is_already_processed=`grep -w '$accession' $tmp_dir\/$continue_from_last_process_file >&1`;
			if($is_already_processed=~/$accession/)
				{
					next;
				}
		}
	#----------------------------------------------------------------------------------------------------


	
	#---------- create a temp file, and start writting the arrays in tab delimited order ----- 
	 my $tmp_output_file=&get_unique_id()."temp_output.txt";
	# system("echo '' >$tmp_dir\/$tmp_output_file");
	 open(WR1,">$tmp_dir\/$tmp_output_file") or print "$!";close(WR1);system("chmod 777 $tmp_dir\/$tmp_output_file");
	 #print "Tmp output in : $tmp_dir\/$tmp_output_file\n\n";
	#----------------------------------------------------------------------------------------- 
	
	#$total_records--;
	 
	 
	 #if($hash_of_already_processed_arrays{$accession}){next;}
	
	 if(not -e "$tmp_dir\/$accession\.fna")
		{
			print "\tNot found \$accession=$accession\n";
			next;
		}
	 
	 #print "\$accession=$accession\n";#next;
	 
	#$pm->start and next; # do the fork 
	 
	 


	#$accession.time().".gff"; 
	 
	#print "Current accession=$accession \tRemaining: $total_records\n\n<br>";	#next;
	#if($accession!~/NC_000913/){exit 0; next;}
	
	my %hash_of_arrays_per_accession;
	my %hash_of_original_arrays_per_accession;
	#my %hash_of_questionable_arrays_per_accession;
	
	#-------------------------------------------------------------------
		
	my $crispr_index=1;
	my $processed_hotspots=0;
	foreach my $range( keys %{$hash_of_all_crispr_hotspots{$accession}})
	  {
		 $total_records--;
		 
		 #my($r_start,$r_stop)=split("-",$range);	
		 #if($r_start<3854077 or $r_start>3854521){next;}
		 
		 
		 #sleep(1);   #----- don't block this, or else ./water will create problem
		 #usleep(1000);
		 select(undef, undef, undef, 0.05); #--- will sleep for 1/4 seconds
		 #-----------------
		 $pm->start and next; # do the fork 
		 
		 
		  
		my $not_a_crispr=0;  		
		my $remaining_hotspots=$total_records-$processed_hotspots;
		
		$processed_hotspots++;
		
		#print "\n\n";  
		
		if($quiet !=1)
			{
				if($hash_id_lookup_table{$accession})
					{
						print "\t$accession\t$hash_id_lookup_table{$accession} - $range\t Remaining: $remaining_hotspots\n";
					}
				else{
						print "\t$accession\t - $range\n";
					}	
			}	
		#---------------
		
		
		  

		
		#------------------------- store the whole array in @current_array after fixing the coords  ---------------------------------------------------------------------------
		my($range_start,$range_stop)=split("-",$range);		
		$hash_of_original_arrays_per_accession{$range_start}{$range_stop}=$hash_of_all_crispr_hotspots{$accession}{$range};
		
		
		
		
		
		
		
		
		
		#----- check if the $range_start already found within another array --------------------------------------
		

		#------ --------------------------------------------------------------------------------------------------
		my ($crispr_hotspot_already_exist,$existing_q_score)=&check_array_existance($array_quality_score_cutoff,$range_start,$range_stop,$tmp_output_file);
		
		
		my $skip_c1=0;
		if($skip_c1==1)
		{
		#foreach my $q_score(sort{$b<=>$a} keys %hash_of_arrays_per_accession)
		#{
		foreach my $asp(sort{$a<=>$b} keys %hash_of_arrays_per_accession)
			{
				#print "\$asp=$asp\n";
				foreach my $astp(sort{$a<=>$b} keys %{$hash_of_arrays_per_accession{$asp}})
					{	
						#print "\$astp=$astp\n";
						
						my $lb; my $ub;
						if($asp<$astp){$lb=$asp;$ub=$astp;}
						else{$lb=$astp;$ub=$asp;}
						
						my $array_center=int($lb+($ub-$lb)/2);
						my $existing_middle_point=int($lb+($ub-$lb)/2);
						my $current_middle_point=int($range_start+($range_stop-$range_start)/2);
						
						#print "START: $asp\tSTOP: $astp\t LB:$lb\tUB:$ub \t\$array_start_position=$array_start_position\t\$array_stop_position=$array_stop_position\n";
						
						#---- check if either $array_start_position or $array_stop_position matches with existing record or not --------------------------------------
						
						if(($current_middle_point>$lb and $current_middle_point<$ub) or ($existing_middle_point>$range_start and $existing_middle_point<$range_stop))
							{
								$crispr_hotspot_already_exist=1;
								last;
							}
					}
				if($crispr_hotspot_already_exist==1){last;}	
			}
			#if($crispr_hotspot_already_exist==1){last;}
		#}					
		}
		
		
		
		if($crispr_hotspot_already_exist==1)
			{
				if($quiet !=1)
					{
						#print "\t\tAlready part of previously processed CRISPR. $range Skipping...\n";
					}	
				#next;
				$pm->finish; 
			}
		
		
		my @original_array=split('\n',$hash_of_all_crispr_hotspots{$accession}{$range});
		
		#-------------- get the model repeat ---------------------------------------------------------	
		my $model_repeat;	
		my $last_line=$original_array[$#original_array];$last_line=~s/\s+/\t/g;
		my @tmp_arr_last_line=split('\t',$last_line);
		$model_repeat=$tmp_arr_last_line[$#tmp_arr_last_line];
		#---------------------------------------------------------------------------------------------
		#foreach my $l1(@original_array)
		#	{
		#		print "$l1\n";
		#	}
		#next;
		
		#my $first_occurrence_of_gap=0;
		#my $coord_diff=0;
		
		
		

		

		
		#--------------------------------------------------- get @current_array ----------------------

		my @current_array; 
		my @modified_array;
		
		my($avg_spacer_length,$total_spacer_length)=&get_current_array(\@original_array,\@current_array);
		

			
		#print "\n\nInput array:\n";		
		#&print_array($model_repeat,\@current_array);
		#---------------------------------------------
		
	my $atleast_one_operation_performed=0;

	#-------------- Now apply different modules -----------------------------------
			
	for(my $pass=0;$pass<2;$pass++)
		{	
			my($array_start_pos,$array_stop_pos,$avg_spacer_len)=&find_array_start_stop_position('F',\@current_array);			
			$avg_spacer_length=$avg_spacer_len;
			
			#print "Pass: $pass $array_start_pos,$array_stop_pos,$avg_spacer_len\n";			
			undef @modified_array;
			
			#----------Special case: remove poor/degenerated repeats from either ends ------------
			
			my $no_of_repeats=$#current_array-4;
			
			if($pass==1)
				{							
					
					
					my $first_repeat_line=$current_array[4];
					my @arr_t1=split('\t',$first_repeat_line);
					my $first_repeat_seq=$arr_t1[1];
					my $tmp_frs=$first_repeat_seq;$tmp_frs=~s/\.//g;
					my %tmp_hash_frs;
					my $number_of_insertions_frs=0;
					if(defined $arr_t1[3]){ $number_of_insertions_frs=&get_number_of_insertions($arr_t1[3],\%tmp_hash_frs);}
					my $first_repeat_degeneracy_score=length($tmp_frs)+$number_of_insertions_frs;
					
					
					
					my $last_repeat_line=$current_array[$#current_array-1];
					my @arr_t2=split('\t',$last_repeat_line);
					my $last_repeat_seq=$arr_t2[1];
					my $tmp_lrs=$last_repeat_seq;$tmp_lrs=~s/\.//g;
					my %tmp_hash_lrs;
					my $number_of_insertions_lrs=0;
					if(defined $arr_t2[3]){$number_of_insertions_lrs=&get_number_of_insertions($arr_t2[3],\%tmp_hash_lrs);}
					my $last_repeat_degeneracy_score=length($tmp_lrs)+$number_of_insertions_lrs;
					
					#print "\$no_of_repeats=$no_of_repeats\t \$first_repeat_seq=$first_repeat_seq\t \$last_repeat_seq=$last_repeat_seq\n";
					
					my @side_to_shorten_first;
					if($no_of_repeats>3 and $first_repeat_degeneracy_score!=$last_repeat_degeneracy_score and ($first_repeat_degeneracy_score>0 or $last_repeat_degeneracy_score>0))
						{
							if($first_repeat_degeneracy_score>=$last_repeat_degeneracy_score )
								{
									push(@side_to_shorten_first,"TOP-1,");
									#push(@side_to_shorten_first,"BOTTOM-1,");
									push(@side_to_shorten_first,"NA,");
								}
							else{
									push(@side_to_shorten_first,"BOTTOM-1,");
									#push(@side_to_shorten_first,"TOP-1,");
									push(@side_to_shorten_first,"NA,");
								}	
						
						
					

					if($shorten_array==1)
							{
									
									foreach my $tmp_user_side_to_shorten(@side_to_shorten_first)
										{
									
											my $case_found=0;	
											my $tmp_minimum_no_of_repeats=1;
											my $tmp_sa_allowed_percent_similarity=90;							

											#------ finally remove the rest degenerated repeats		
											#$tmp_user_side_to_shorten="BOTTOM-1,";
											($case_found)=&shorten_array($range,$tmp_minimum_no_of_repeats,$tmp_sa_allowed_percent_similarity,$tmp_user_side_to_shorten,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
									
											if($case_found==1)
												{
													for(my $m=3;$m>=0;$m--)
														{
															unshift(@modified_array,$current_array[$m]);
															#shift(@modified_array);
														}
													push(@modified_array,$current_array[$#current_array]);	
													
													@current_array=@modified_array;
													$atleast_one_operation_performed++;
												}
											
											
											undef @modified_array;
										}	

							}
							
					if($remove_insertion_from_repeat==1 and $model_repeat=~/-/)
							{

									
									my $case_found=0;
									($model_repeat,$case_found)=&fix_arrays_with_insertion_in_repeat($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
									
									
									if($case_found==1)
										{
											$atleast_one_operation_performed++;
										}
									


									undef @modified_array;	
							}
							
							
					if($check_consensus==1)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
							{
									#print "\nA1: $left_flank_seq\nA2: $right_flank_seq\n\n";
									
									my $case_found=0;	
									($model_repeat,$case_found)=&check_consensus_sequence($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
									
									
									if($case_found==1)
										{
											$atleast_one_operation_performed++;
										}
									undef @modified_array;				
							}
									

							#print "$pass:Array after shorten_array\n";
							#&print_array($model_repeat,\@current_array);	
						}							
				}




			#--------------------------------- this is default flow [the first two modules has the be in the order: 1. Extend repeat length, then trim repeat length]--------------------------------------------------------------------------

			#if($pass==0)
			#	{
				
					######################### [DO not change the flow and order of the first 2 modules: 1increase_repeat_length, 2: trim_repeats, 3: extend
					
					#print "A1: $pass:Before increase_repeat_length\n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------			
					
					if($increase_repeat_length==1)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
						{
							#system("echo '$user_side_to_increase_length' >log.txt");
							my @arr_user_side_to_increase_length=split(',',$user_side_to_increase_length);
							
							foreach my $side_and_bases(@arr_user_side_to_increase_length)
								{			
							
									if($side_and_bases!~/NA/)
										{
											my($u_side_to_increase_length,$u_no_of_bases_to_increase_length)=split('-',$side_and_bases);
											#---------------------------------------------------------------
											
											my $case_found=0;	
											undef @modified_array;
											($model_repeat,$case_found)=&increase_all_repeat_lengths($range,$u_side_to_increase_length,$u_no_of_bases_to_increase_length,$accession,$model_repeat,\@current_array,\@modified_array);
											
											
											if($case_found==1)
												{
													#------------ complete the @modified array ---
													for(my $m=3;$m>=0;$m--)
														{
															unshift(@modified_array,$current_array[$m]);
															#shift(@modified_array);
														}
													push(@modified_array,$current_array[$#current_array]);
														
													@current_array=@modified_array;
													
													#undef @modified_array;
													$atleast_one_operation_performed++;
												}
												
											undef @modified_array;	
										}	
								
									else{							#----- auto check with $repeat_extension_identity
									
											#system("echo '$repeat_extension_identity' >log.txt");
											
											my $case_found1=0;	
											my $sides_to_increase_length="NA-0,";
											undef @modified_array;
											
											($sides_to_increase_length,$case_found1)=&auto_detect_repeat_sides_and_bases_to_increase($range,$repeat_extension_identity,$accession,$model_repeat,\@current_array,\@modified_array);
											
											#print "$sides_to_increase_length,$case_found1\n";
											
											if($case_found1==1)
												{
													
													my @arr_sides_to_increase_length=split(',',$sides_to_increase_length);									
													
													foreach my $side_and_bases(@arr_sides_to_increase_length)
														{			
															#system("echo '$side_and_bases' >log.txt");
															
															if($side_and_bases!~/NA/)
																{
																	my($side_to_increase_length,$no_of_bases_to_increase_length)=split('-',$side_and_bases);
																	#---------------------------------------------------------------
																	
																	my $case_found=0;	
																	#print "$range,$side_to_increase_length,$no_of_bases_to_increase_length,$accession,$model_repeat,\n";
																	($model_repeat,$case_found)=&increase_all_repeat_lengths($range,$side_to_increase_length,$no_of_bases_to_increase_length,$accession,$model_repeat,\@current_array,\@modified_array);
																	
																	
																	if($case_found==1)
																		{
																			#------------ complete the @modified array ---
																			for(my $m=3;$m>=0;$m--)
																				{
																					unshift(@modified_array,$current_array[$m]);
																					#shift(@modified_array);
																				}
																			push(@modified_array,$current_array[$#current_array]);
																				
																			@current_array=@modified_array;
																			
																			#undef @modified_array;
																			$atleast_one_operation_performed++;
																		}
																		
																	undef @modified_array;	
																}
														}		
												}
												
											undef @modified_array;	
										}	
								
								}	
						}
					
					
					#print "A2: $pass:After increase_repeat_length\n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------
					
					
					if($trim_repeats==1)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
						{							
							my @arr_user_sides_to_trim=split(',',$user_side_to_trim);
							
							foreach my $side_and_bases(@arr_user_sides_to_trim)
								{
									
									my $case_found=0;				
									#$minimum_repeat_length=11;
									my $trimming_cutoff2=$trimming_cutoff;
									
									if($no_of_repeats>3){$trimming_cutoff2=20;}
									
									($model_repeat,$case_found)=&trim_repeat_ends_and_improve_score($range,$minimum_repeat_length,$side_and_bases,$trimming_cutoff2,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
									
									
									if($case_found==1)
										{
											$atleast_one_operation_performed++;
										}
									
									
									if($case_found==1)
										{	
											
											if($side_and_bases=~/NA/)
												{
													undef @modified_array;							
													for(my $m=3;$m>=0;$m--)
														{
															unshift(@modified_array,$current_array[$m]);
														}
													
													my $case_found1=0;
													($case_found1)=&fix_arrays_with_gaps($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
													
													
													push(@modified_array,$current_array[$#current_array]);
																					
													
													if($case_found1==1)	{@current_array=@modified_array;}								
																
													undef @modified_array;	
												}
										
									}
										
											
									undef @modified_array;
								}
							
						}
					
					#print "A2.1: $pass:after trim_repeats\n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------
					
					if(length($model_repeat)<=$word_length){$not_a_crispr=1; last;}
				
										
					
					
					
					#&fix_array_with_falsely_identified_repeats($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);	
					
					#---------------------------------------------------------
					#print "A3: $pass:After fix_array_with_falsely_identified_repeats\n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------
					
					
					#&check_array_for_misaligned_bases($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
					
					#---------------------------------------------------------
					#print "A4: $pass:After check_array_for_misaligned_bases \n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------
					
					if($remove_insertion_from_repeat==1 and $model_repeat=~/-/)
						{
							
							my $case_found=0;
							($model_repeat,$case_found)=&fix_arrays_with_insertion_in_repeat($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
							
							
							if($case_found==1)
								{
									$atleast_one_operation_performed++;
								}
							
	
							undef @modified_array;	
						}
					#next;
					
				
					if($check_consensus==1)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
						{
							#print "\nA1: $left_flank_seq\nA2: $right_flank_seq\n\n";
							
							my $case_found=0;	
							($model_repeat,$case_found)=&check_consensus_sequence($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
							
							
							if($case_found==1)
								{
									$atleast_one_operation_performed++;
								}
							undef @modified_array;				
						}
						
					#---------------------------------------------------------
					#print "A5: $pass:After check_consensus \n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------
					

				
				
					
					#---------- try to extend arrays in flanks -----------------------------------
					
					#print "Going to check extension\n";
					
					if($extend_array==1 and length($model_repeat)>=$word_length)    #----- remember; don't use ea_dynamic search at the first extension
						{
							
							my $case_found=0;	
							#--- for the first extension dont use ea_dynamic_search
							my $tmp_ea_dynamic_search=0;		
							($case_found)=&extend_array($range,$max_gap_between_crisprs,$tmp_ea_dynamic_search,$ea_allowed_percent_similarity,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
							
							
							if($case_found==1)
								{
									$atleast_one_operation_performed++;
								}
							
							if($case_found==1)
								{
									for(my $m=3;$m>=0;$m--)
										{
											unshift(@modified_array,$current_array[$m]);
											#shift(@modified_array);
										}
									push(@modified_array,$current_array[$#current_array]);	
											
									@current_array=@modified_array;
									undef @modified_array;
									
									
									#------- fix the gaps ---------------------------------
									for(my $m=3;$m>=0;$m--)
										{
											unshift(@modified_array,$current_array[$m]);
										}
									
									my $case_found1=0;
									($case_found1)=&fix_arrays_with_gaps($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
									
								
									if($case_found1==1)
										{
											push(@modified_array,$current_array[$#current_array]);
											@current_array=@modified_array;	
											$atleast_one_operation_performed++;
										}
												
									undef @modified_array;	
									
								}
							
							
							undef @modified_array;
						}
					

					#print "A7: $pass:After extend_array\n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------	

					if($remove_insertion_from_repeat==1 and $model_repeat=~/-/)
						{

							
							my $case_found=0;
							($model_repeat,$case_found)=&fix_arrays_with_insertion_in_repeat($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
							
							
							if($case_found==1)
								{
									$atleast_one_operation_performed++;
								}
							
	
							undef @modified_array;	
						}
					#next;
					
				
					if($check_consensus==10)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
						{
							#print "\nA1: $left_flank_seq\nA2: $right_flank_seq\n\n";
							
							my $case_found=0;	
							($model_repeat,$case_found)=&check_consensus_sequence($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
							
							
							if($case_found==1)
								{
									$atleast_one_operation_performed++;
								}
							undef @modified_array;				
						}
					
				
					#print "A8: $pass:After check_consensus\n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------				
					

					if(($#current_array-4)>3)
						{												
							if($trim_repeats==1)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
								{									
									my @arr_user_sides_to_trim=split(',',$user_side_to_trim);
									
									foreach my $side_and_bases(@arr_user_sides_to_trim)
										{
											
											my $case_found=0;				
											#$minimum_repeat_length=11;
											($model_repeat,$case_found)=&trim_repeat_ends_and_improve_score($range,$minimum_repeat_length,$side_and_bases,$trimming_cutoff,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
											
											#print "\t$model_repeat,$case_found\n";
											if($case_found==1)
												{
													$atleast_one_operation_performed++;
												}
											
											
											if($case_found==1)
												{	
													
													if($side_and_bases=~/NA/)
														{
															undef @modified_array;							
															for(my $m=3;$m>=0;$m--)
																{
																	unshift(@modified_array,$current_array[$m]);
																}
															
															my $case_found1=0;
															($case_found1)=&fix_arrays_with_gaps($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
															
															
															push(@modified_array,$current_array[$#current_array]);
																							
															
															if($case_found1==1)	{@current_array=@modified_array;}								
																		
															undef @modified_array;	
														}

											}
												
													
											undef @modified_array;
										}									
								}	
								
							#print "A8.1: $pass:After ($#current_array-4)>3\t trim_repeats\n";
							#&print_array($model_repeat,\@current_array);
							#----------------------------------------------------------		
							
							#------------------------------------------------------------
							
							if($remove_insertion_from_repeat==1 and $model_repeat=~/-/)
								{

									
									my $case_found=0;
									($model_repeat,$case_found)=&fix_arrays_with_insertion_in_repeat($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
									
									
									if($case_found==1)
										{
											$atleast_one_operation_performed++;
										}
									

									undef @modified_array;	
								}
							
							
							if($check_consensus==1)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
								{
									#print "\nA1: $left_flank_seq\nA2: $right_flank_seq\n\n";
									
									my $case_found=0;	
									($model_repeat,$case_found)=&check_consensus_sequence($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
									
									
									if($case_found==1)
										{
											$atleast_one_operation_performed++;
										}
									undef @modified_array;				
								}
			
						}			

					#print "A9: $pass:After ($#current_array-4)>3\n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------	
					
					#---------- then remove poor/degenerated repeats from either ends ------------

					if($shorten_array==1)
						{
							

							my $case_found=0;	
							#my $tmp_minimum_no_of_repeats=2;
							#my $tmp_sa_allowed_percent_similarity=80;		
							($case_found)=&shorten_array($range,$minimum_no_of_repeats,$sa_allowed_percent_similarity,$user_side_to_shorten,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
							
							

							
							if($case_found==1)
								{
									for(my $m=3;$m>=0;$m--)
										{
											unshift(@modified_array,$current_array[$m]);
											#shift(@modified_array);
										}
									push(@modified_array,$current_array[$#current_array]);	
									
									
									#print "\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n";	
									#for(my $k1=0;$k1<=$#original_array;$k1++)
									#	{
									#		print "$original_array[$k1]\n";
									#	} 			
									#print "\n";
									
									#print "\n1.LF: $left_flank_seq\nRF: $right_flank_seq\n";
									
									#print "Searched for extension of the array in flanks of $accession:\n\n";
											#print "\@modified_array=@modified_array\n";
									
									
									#print "Position\tRepeat\tSpacer\tInsertion\n";
									#print "========\t======\t======\t=========\n";
									
									#foreach my $l(@modified_array)
									#	{
									#		print "A: $l\n";
									#	}
									#print "\$model_repeat=$model_repeat\n\n";	
									
									
									#print "\n";
									
									@current_array=@modified_array;
									$atleast_one_operation_performed++;
								}
							
							
							undef @modified_array;
						}

					#print "A10: $pass:After shorten_array\n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------						
					
					if($remove_insertion_from_repeat==1 and $model_repeat=~/-/)
						{							
							my $case_found=0;
							($model_repeat,$case_found)=&fix_arrays_with_insertion_in_repeat($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
							
							
							if($case_found==1)
								{
									$atleast_one_operation_performed++;
								}
							
	
							undef @modified_array;	
						}
					
					#print "$pass:After remove_insertion_from_repeat\n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------		
					
					if($check_consensus==1)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
						{
							#print "\nA1: $left_flank_seq\nA2: $right_flank_seq\n\n";
							
							my $case_found=0;	
							($model_repeat,$case_found)=&check_consensus_sequence($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
							
							
							if($case_found==1)
								{
									$atleast_one_operation_performed++;
								}
							undef @modified_array;				
						}
							

					#print "A11: $pass:After check_consensus\n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------					


					
					if($trim_repeats==1)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
						{
							
							
							my @arr_user_sides_to_trim=split(',',$user_side_to_trim);
							
							foreach my $side_and_bases(@arr_user_sides_to_trim)
								{
									
									my $case_found=0;				
									#$minimum_repeat_length=11;
									($model_repeat,$case_found)=&trim_repeat_ends_and_improve_score($range,$minimum_repeat_length,$side_and_bases,$trimming_cutoff,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
									
									
									if($case_found==1)
										{
											$atleast_one_operation_performed++;
										}
									
									
									if($case_found==1)
										{	
											
											if($side_and_bases=~/NA/)
												{
													undef @modified_array;							
													for(my $m=3;$m>=0;$m--)
														{
															unshift(@modified_array,$current_array[$m]);
														}
													
													my $case_found1=0;
													($case_found1)=&fix_arrays_with_gaps($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
													
													
													push(@modified_array,$current_array[$#current_array]);
																					
													
													if($case_found1==1)	{@current_array=@modified_array;}								
																
													undef @modified_array;	
												}
											#print "\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n";	
											#for(my $k1=0;$k1<=$#original_array;$k1++)
											#	{
											#		print "$original_array[$k1]\n";
											#	} 			
											#print "\n";
											
											#print "Repeat ends trimming operation for $accession:\n\n";
											
											#foreach my $l(@current_array)
											#	{
											#		print "$l\n";
											#	}
											#print "\$model_repeat=$model_repeat\n\n";
											
											#print "\nLF: $left_flank_seq\nRF: $right_flank_seq\n";	
											
											#--- check if the insertion is already handled and removed or not
											#if($model_repeat!~/-/)
											#	{
											#		#$insertion_found_in_repeat=0;
											#	}
									}
										
											
									undef @modified_array;
								}
							
						}
					
					#print "A12: $pass:After trim_repeats\n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------		
					
					
					if($remove_insertion_from_repeat==1 and $model_repeat=~/-/)
						{						
							my $case_found=0;
							($model_repeat,$case_found)=&fix_arrays_with_insertion_in_repeat($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
							
							
							if($case_found==1)
								{
									$atleast_one_operation_performed++;
								}

							undef @modified_array;	
						}
					
					
					if($check_consensus==1)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
						{
							#print "\nA1: $left_flank_seq\nA2: $right_flank_seq\n\n";
							
							my $case_found=0;	
							($model_repeat,$case_found)=&check_consensus_sequence($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
							
							
							if($case_found==1)
								{
									$atleast_one_operation_performed++;
								}
							undef @modified_array;				
						}
					
					#print "A13: $pass:After trim_repeats\n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------				
			#	}
		



					#print "\n\n\nAAAAAAAAA:::::: $pass:Array after pass:$pass\n";
					#&print_array($model_repeat,\@current_array);
					#----------------------------------------------------------	


					########################### if the median spacer length become shorter than 5, give up ############################	
					
					my $median_spacer_length=&get_median_spacer_length(\@current_array);	
					
					#print "\$median_spacer_length=$median_spacer_length\n";
					if($median_spacer_length<5)
						{
							$not_a_crispr=1;
							#$pm->finish;							
							#last;						
						}
					################################################################################################################

	} #------- end of pass		
		
		

#----- skip to next record, if $model_repeat_length <15

		###################################################################################################################
		#----------------------------------------- post processing modules -----------------------------------------------#
		###################################################################################################################		
		
		
		
	my $median_spacer_length=&get_median_spacer_length(\@current_array);

	
	if($not_a_crispr==0 and length($model_repeat)>$repeat_length_cutoff and $median_spacer_length>17)
		{		
			#print "Before search_unidentified_repeat_in_spacer_sequence\n";
			#&print_array($model_repeat,\@current_array);
			#----------------------------------------------------------		
			
			#my $no_of_repeats_in_array=$#current_array-4;
			#print "$#current_array-4\n";
			#---------- try to extend arrays in flanks -----------------------------------
					
					
			if($extend_array==1 and ($#current_array-4) < 4 and length($model_repeat)>=$word_length)
				{
								
								my $case_found=0;	
										
								($case_found)=&extend_array($range,$max_gap_between_crisprs,$ea_dynamic_search,$ea_allowed_percent_similarity,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
								
								
								if($case_found==1)
									{
										$atleast_one_operation_performed++;
									}
								
								if($case_found==1)
									{
										for(my $m=3;$m>=0;$m--)
											{
												unshift(@modified_array,$current_array[$m]);
												#shift(@modified_array);
											}
										push(@modified_array,$current_array[$#current_array]);	
												
										@current_array=@modified_array;
										undef @modified_array;
										
										
										#------- fix the gaps ---------------------------------
										for(my $m=3;$m>=0;$m--)
											{
												unshift(@modified_array,$current_array[$m]);
											}
										
										my $case_found1=0;
										($case_found1)=&fix_arrays_with_gaps($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
										
									
										if($case_found1==1)
											{
												push(@modified_array,$current_array[$#current_array]);
												@current_array=@modified_array;	
												$atleast_one_operation_performed++;
											}
													
										undef @modified_array;	
										
									}
								
								
								undef @modified_array;
				}
						

		
			if($check_consensus==1)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
				{
					#print "\nA1: $left_flank_seq\nA2: $right_flank_seq\n\n";
					
					my $case_found=0;	
					($model_repeat,$case_found)=&check_consensus_sequence($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
					
					
					if($case_found==1)
						{
							$atleast_one_operation_performed++;
						}
					undef @modified_array;				
				}
		
		
			#print "Before search_unidentified_repeat_in_spacer_sequence\n";
			#&print_array($model_repeat,\@current_array);
			#----------------------------------------------------------	
					
				#print "Going to check search_unidentified_repeat_in_spacer_sequence\n";		
			if($search_unidentified_repeat_in_spacer_sequence==1)
				{
					
					my($array_start_pos,$array_stop_pos,$avg_spacer_len)=&find_array_start_stop_position('F',\@current_array);			
					$avg_spacer_length=$avg_spacer_len;
					#print "$array_start_pos,$array_stop_pos,$avg_spacer_len\n\n";
					
					for(my $m=3;$m>=0;$m--)
						{
							unshift(@modified_array,$current_array[$m]);
						}
					my $case_found=0;
					#my $allowed_percent_similarity=60;				
					($case_found)=&search_unidentified_repeat_in_spacers($range,$su_dynamic_search,$allowed_percent_similarity,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);				
					
					push(@modified_array,$current_array[$#current_array]);
					
					
					
					if($case_found==1)
						{
							$atleast_one_operation_performed++;
						}
						
					if($case_found==1)
						{
							@current_array=@modified_array;
						}
					
					
					undef @modified_array;
				}
				
			#print "After search_unidentified_repeat_in_spacer_sequence\n";
			#&print_array($model_repeat,\@current_array);
		
		


		

		if($fix_gaps_in_repeats==1)   # this will be called repeat_end_correction module
			{
				
				
				my $case_found=0;
				($case_found)=&fix_arrays_with_gaps($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
				
					
				
				if($case_found==1)
					{
						for(my $m=3;$m>=0;$m--)
							{
								unshift(@modified_array,$current_array[$m]);
							}
						
						push(@modified_array,$current_array[$#current_array]);
						@current_array=@modified_array;	
						$atleast_one_operation_performed++;
					}
							
				undef @modified_array;	
					
			}		




		if($remove_insertion_from_repeat==1 and $model_repeat=~/-/)
			{

				
				my $case_found=0;
				($model_repeat,$case_found)=&fix_arrays_with_insertion_in_repeat($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
				
				
				if($case_found==1)
					{
						$atleast_one_operation_performed++;
					}				

				if($case_found==1)
					{

					}	
				undef @modified_array;	
			}
		#next;
		
		if($check_consensus==1)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
			{
				#print "\nA1: $left_flank_seq\nA2: $right_flank_seq\n\n";
				
				my $case_found=0;	
				($model_repeat,$case_found)=&check_consensus_sequence($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
				
				
				if($case_found==1)
					{
						$atleast_one_operation_performed++;
					}
				undef @modified_array;				
			}




		########################################## special cases: $repeat_extension_identity set to 75 (and $alt_repeat_extension_identity=50) ####

		if($increase_repeat_length==1 and length($model_repeat)<24)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
			{
				#system("echo '$user_side_to_increase_length' >log.txt");
				my @arr_user_side_to_increase_length=split(',',$user_side_to_increase_length);
				
				foreach my $side_and_bases(@arr_user_side_to_increase_length)
					{			
				
						if($side_and_bases!~/NA/)
							{
								my($u_side_to_increase_length,$u_no_of_bases_to_increase_length)=split('-',$side_and_bases);
								#---------------------------------------------------------------
								
								my $case_found=0;	
								undef @modified_array;
								($model_repeat,$case_found)=&increase_all_repeat_lengths($range,$u_side_to_increase_length,$u_no_of_bases_to_increase_length,$accession,$model_repeat,\@current_array,\@modified_array);
								
								
								if($case_found==1)
									{
										#------------ complete the @modified array ---
										for(my $m=3;$m>=0;$m--)
											{
												unshift(@modified_array,$current_array[$m]);
												#shift(@modified_array);
											}
										push(@modified_array,$current_array[$#current_array]);
											
										@current_array=@modified_array;
										
										#undef @modified_array;
										$atleast_one_operation_performed++;
									}
									
								undef @modified_array;	
							}	
					
						else{							#----- auto check with $repeat_extension_identity
						
								#system("echo '$repeat_extension_identity' >log.txt");
								
								my $case_found1=0;	
								my $sides_to_increase_length="NA-0,";
								undef @modified_array;
								
								($sides_to_increase_length,$case_found1)=&auto_detect_repeat_sides_and_bases_to_increase($range,$repeat_extension_identity,$accession,$model_repeat,\@current_array,\@modified_array); # special cases: $repeat_extension_identity set to 50
								
								#print "$sides_to_increase_length,$case_found1\n";
								
								if($case_found1==1)
									{
										
										my @arr_sides_to_increase_length=split(',',$sides_to_increase_length);									
										
										foreach my $side_and_bases(@arr_sides_to_increase_length)
											{			
												#system("echo '$side_and_bases' >log.txt");
												
												if($side_and_bases!~/NA/)
													{
														my($side_to_increase_length,$no_of_bases_to_increase_length)=split('-',$side_and_bases);
														#---------------------------------------------------------------
														
														my $case_found=0;	
														
														($model_repeat,$case_found)=&increase_all_repeat_lengths($range,$side_to_increase_length,$no_of_bases_to_increase_length,$accession,$model_repeat,\@current_array,\@modified_array);
														
														
														if($case_found==1)
															{
																#------------ complete the @modified array ---
																for(my $m=3;$m>=0;$m--)
																	{
																		unshift(@modified_array,$current_array[$m]);
																		#shift(@modified_array);
																	}
																push(@modified_array,$current_array[$#current_array]);
																	
																@current_array=@modified_array;
																
																#undef @modified_array;
																$atleast_one_operation_performed++;
															}
															
														undef @modified_array;	
													}
											}		
									}
									
								undef @modified_array;	
							}	
					
					}	
			
				
				#----- check for misaligned bases
				#&fix_ array_with_falsely_identified_repeats($accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);	 --- use in general
				&check_array_for_misaligned_bases($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
				
			}
		

		if($remove_insertion_from_repeat==1 and $model_repeat=~/-/)
			{
				#print "Going to fix Insertion in repeat:\n\n";
				#print "\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n";
				#for(my $k1=0;$k1<=$#original_array;$k1++)
				#	{
				#		print "$original_array[$k1]\n";
				#	} 
					
				#--- now fix the gaps and print them as we go
				# print "\n\nFixed Array:\n\n";	
				
				#print "A1: $left_flank_seq\nA2: $right_flank_seq\n";
				
				my $case_found=0;
				($model_repeat,$case_found)=&fix_arrays_with_insertion_in_repeat($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
				
				
				if($case_found==1)
					{
						$atleast_one_operation_performed++;
					}
				
				#print "\n";
				
						#print "\@modified_array=@modified_array\n";
				#print "Position\tRepeat\tSpacer\tInsertion\n";
				#print "========\t======\t======\t=========\n";
				if($case_found==1)
					{
						
						#print "\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n";	
						#for(my $k1=0;$k1<=$#original_array;$k1++)
						#	{
						#		print "$original_array[$k1]\n";
						#	}
						
						#print "Searched for insertion in model_repeat of $accession:\n\n";
						
						#foreach my $l(@current_array)
						#	{
						#		print "$l\n";
						#	}
						#print "\nLF: $left_flank_seq\nRF: $right_flank_seq\n";	
						#print "\$model_repeat=$model_repeat\n\n";
						
						#$insertion_found_in_repeat=0;
					}	
				undef @modified_array;	
			}
		#next;
		
		if($check_consensus==1)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
			{
				#print "\nA1: $left_flank_seq\nA2: $right_flank_seq\n\n";
				
				my $case_found=0;	
				($model_repeat,$case_found)=&check_consensus_sequence($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
				
				
				if($case_found==1)
					{
						$atleast_one_operation_performed++;
					}
				undef @modified_array;				
			}
		


		if($trim_repeats==1 and length($model_repeat)>24)######### $trimming_cutoff set to high ( 20 ) ; $minimum_repeat_length set to 24
			{
				
				
				my @arr_user_sides_to_trim=split(',',$user_side_to_trim);
				
				foreach my $side_and_bases(@arr_user_sides_to_trim)
					{
						
						my $case_found=0;				
						my $tmp_minimum_repeat_length=24;
						($model_repeat,$case_found)=&trim_repeat_ends_and_improve_score($range,$tmp_minimum_repeat_length,$side_and_bases,$trimming_cutoff,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
						
						
						if($case_found==1)
							{
								$atleast_one_operation_performed++;
							}
						
						
						if($case_found==1)
							{	
								
								if($side_and_bases=~/NA/)
									{
										undef @modified_array;							
										for(my $m=3;$m>=0;$m--)
											{
												unshift(@modified_array,$current_array[$m]);
											}
										
										my $case_found1=0;
										($case_found1)=&fix_arrays_with_gaps($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
										
										
										push(@modified_array,$current_array[$#current_array]);
																		
										
										if($case_found1==1)	{@current_array=@modified_array;}								
													
										undef @modified_array;	
									}
								#print "\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n";	
								#for(my $k1=0;$k1<=$#original_array;$k1++)
								#	{
								#		print "$original_array[$k1]\n";
								#	} 			
								#print "\n";
								
								#print "Repeat ends trimming operation for $accession:\n\n";
								
								#foreach my $l(@current_array)
								#	{
								#		print "$l\n";
								#	}
								#print "\$model_repeat=$model_repeat\n\n";
								
								#print "\nLF: $left_flank_seq\nRF: $right_flank_seq\n";	
								
								#--- check if the insertion is already handled and removed or not
								#if($model_repeat!~/-/)
								#	{
								#		#$insertion_found_in_repeat=0;
								#	}
						}
							
								
						undef @modified_array;
					}
				
			}
	
		#print "2: After everything:\n";		
		#&print_array($model_repeat,\@current_array);
		#----------------------------------------------------------	
	
		###################################### END of Special cases ################################################################################




		&fix_array_with_falsely_identified_repeats($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);	

		#----------- check potential alternate repeat -----------------------------------------	
		#print "Going to check potential alternate repeat\n";
			
		if($find_alternate_repeat==1)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
			{
				#print "\nA1: $left_flank_seq\nA2: $right_flank_seq\n\n";
				
				my $case_found=0;	
				($model_repeat,$potential_alternate_repeat,$case_found)=&search_alternate_repeat_sequence($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
				
				
				if($case_found==1)
					{
						$atleast_one_operation_performed++;
					}
					
				undef @modified_array;	
			}

	

		#---------- check direction of the array ---------------------------------------------

		my $direction_found=0;
		my $array_direction="NA";
		my $repeat_family="NA";
		my $matching_reference_repeat="NA";
		my $array_direction_MEMO="NA";
		
		#print "Going to check direction\n";
		
		if($check_direction==1)
			{
							
				
				
				#---- note: array start and stop position is require for Longer leader analysis: find them here which is different from CRISPRDetect.pl --------
				
				my($array_start_position,$array_stop_position,$avg_spacer_length)=&find_array_start_stop_position($array_direction,\@current_array);

				
				#-------------------------------------------------------------------------------
				my $case_found=0;	
				my $matching_reference_repeat_direction="NA";
				
				#($matching_reference_repeat,$model_repeat,$array_direction,$repeat_family,$array_direction_MEMO,$case_found)=&check_array_direction($accession,$model_repeat,$avg_spacer_length,$all_gene_positions_folder,$all_gene_positions_file,\@current_array,\@modified_array,\%lib_of_repeats_with_confirmed_direction);
				($matching_reference_repeat,$matching_reference_repeat_direction,$model_repeat,$array_direction,$repeat_family,$array_direction_MEMO,$case_found)=&check_array_direction($range,$blast_db_file_of_known_repeats,$check_motif_in_repeat,$motif,$check_A_and_T_ratio_in_repeat,$check_similarity_with_reference_repeat,$allowed_no_of_mismatches,$check_secondary_structure_of_repeat,$MFE_cutoff,$MFE_minimum_difference,$MFE_exclude_bases,$check_array_degeneracy,$permitted_mutation_per_array,$check_AT_distribution_in_flanks,$AT_distribution_window,$AT_distribution_minimum_percentage_difference,$check_longer_leader,$Motif_match_score,$A_and_T_ratio_score,$Similarity_score,$MFE_score,$array_degeneracy_score,$AT_distribution_score,$Longer_leader_score,$array_start_position,$array_stop_position,$accession,$model_repeat,$all_gene_positions_folder,$all_gene_positions_file,\%lib_of_repeats_with_confirmed_direction,\@current_array,\@modified_array);
				
				
				my $model_repeat_rc=$model_repeat;$model_repeat_rc=reverse $model_repeat_rc;$model_repeat_rc=~tr/ACGTU/TGCAA/;
				
				if($matching_reference_repeat !~ /^NA/ and ($matching_reference_repeat_direction =~/^$array_direction/ or $array_direction =~/^$matching_reference_repeat_direction/))
					{
						
						#-------- align the model repeat and the reference repeat first, and get the positions of match ---
						my $b_sequence=$matching_reference_repeat;
						if($array_direction=~/R/)
							{
								$b_sequence= reverse $b_sequence;
								$b_sequence=~tr/ACGT/TGCA/;
							}
								
						my($m_repeat_line,$r_repeat_line)=&get_aligned_region_start_stop($range,$accession,$model_repeat,$b_sequence);
								
						#system("echo '$m_repeat_line,$r_repeat_line' >>log1.txt");
						
						
						#---- check $m_repeat_line: to extend the repeats
						if($m_repeat_line=~/^-/ or $m_repeat_line=~/-$/)
							{
								my $side_to_increase_length;
								my $no_of_bases_to_increase_length;
								#-----Check left -----------------------------------------------
								if($m_repeat_line=~/^(\-+)/)
									{									
										$no_of_bases_to_increase_length=length($1);
										$side_to_increase_length="LEFT";
										
										my $case_found1=0;	
										undef @modified_array;
										($model_repeat,$case_found1)=&increase_all_repeat_lengths($range,$side_to_increase_length,$no_of_bases_to_increase_length,$accession,$model_repeat,\@current_array,\@modified_array);
																			
										#------------ complete the @modified array ---
										for(my $m=3;$m>=0;$m--){unshift(@modified_array,$current_array[$m]);}
										push(@modified_array,$current_array[$#current_array]);
													
										@current_array=@modified_array;
									}
								
								#------ check right ----------------------------------------------
								if($m_repeat_line=~/(\-+)$/)
									{									
										$no_of_bases_to_increase_length=length($1);
										$side_to_increase_length="RIGHT";
										
										my $case_found1=0;	
										undef @modified_array;
										($model_repeat,$case_found1)=&increase_all_repeat_lengths($range,$side_to_increase_length,$no_of_bases_to_increase_length,$accession,$model_repeat,\@current_array,\@modified_array);
																			
										#------------ complete the @modified array ---
										for(my $m=3;$m>=0;$m--){unshift(@modified_array,$current_array[$m]);}
										push(@modified_array,$current_array[$#current_array]);
													
										@current_array=@modified_array;
									}
									
								if($model_repeat=~/-/)
									{
										my $case_found_t1=0;	
										($model_repeat,$case_found_t1)=&check_consensus_sequence($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
									}										
								undef @modified_array;	
							}		
						#---- check $r_repeat_line: to trim the repeats
						if($r_repeat_line=~/^-/ or $r_repeat_line=~/-$/)
							{
								my $side_and_bases;										
								my $no_of_bases_to_trim;
										
								#---- first left side ----------------------------------------
								if($r_repeat_line=~/^(\-+)/)
									{
										$no_of_bases_to_trim=length($1);	
										
										
										my $case_found1=0;
										$side_and_bases="LEFT-$no_of_bases_to_trim";
										
										($model_repeat,$case_found1)=&trim_repeat_ends_and_improve_score($range,$minimum_repeat_length,$side_and_bases,$trimming_cutoff,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
												
										undef @modified_array;
									}
								#---- right side ----------------------------------------
								if($r_repeat_line=~/(\-+)$/)
									{
										$no_of_bases_to_trim=length($1);	
										
										
										my $case_found1=0;
										$side_and_bases="RIGHT-$no_of_bases_to_trim";
										
										($model_repeat,$case_found1)=&trim_repeat_ends_and_improve_score($range,$minimum_repeat_length,$side_and_bases,$trimming_cutoff,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
												
										undef @modified_array;
									}	
							}		

												
						$direction_found++;
						$atleast_one_operation_performed++;
					}
				
				elsif($model_repeat=~/ATTGAAA.?\w{0,3}/ or $model_repeat_rc=~/ATTGAAA.?\w{0,3}/)
					{
							my $side_and_bases="NA-0,";										
							my $no_of_bases_to_trim;
										
					
							if($model_repeat=~/ATTGAAA.?(\w{0,3})$/)
								{
									#print "Inside..\n";
									$no_of_bases_to_trim=length($1);
									$side_and_bases="RIGHT-$no_of_bases_to_trim";	
								}		
							elsif($model_repeat_rc=~/ATTGAAA.?(\w{0,3})$/)
								{
									$no_of_bases_to_trim=length($1)-1; # left side starts from zero
									$side_and_bases="RIGHT-$no_of_bases_to_trim";
								}			
							
							
							
							if($side_and_bases!~/NA/)
								{
									my $case_found1=0;										
										
									($model_repeat,$case_found1)=&trim_repeat_ends_and_improve_score($range,$minimum_repeat_length,$side_and_bases,$trimming_cutoff,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
												
									undef @modified_array;
								}	
								
								
					}
				
				
				
				
				
				
				
				
				#--- if $array_direction="R" then reverse complement all the repeats and spacers	
				#if($array_direction =~/R/)
				#	{
				#		#---- if $potential_alternate_repeat exists, then reverse-comp that too
				#		
				#		if(defined $potential_alternate_repeat and $potential_alternate_repeat ne "NA")
				#			{
				#				$potential_alternate_repeat=reverse $potential_alternate_repeat;
				#				$potential_alternate_repeat=~tr/ACGT/TGCA/;
				#			}
				#	}	
				undef @modified_array;	
			}

		
		&check_array_for_misaligned_bases($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
		#print "Done\n";
		



		#---------- get array_quality score
		#my $species="Bacteria";
		if($hash_of_accession_and_species{$accession})
			{
				$species=$hash_of_accession_and_species{$accession};
			}
		my($array_quality_score,$score_det,$score_legend)=&calculate_array_quality_score($range,$species,$accession,$all_gene_positions_folder,$all_gene_positions_file,$matching_reference_repeat,$model_repeat,\@current_array);
		$array_quality_score= sprintf("%.2f",$array_quality_score);
		
	
		
		#print "$range: \$array_quality_score=$array_quality_score\n";
		#------- double check the arrays, where score is close to 0
		if($array_quality_score<0 and $array_quality_score > -1.50)   ### remove degenerated repeats from either end to improve score
			{
				#---------- remove poor/degenerated repeats from either ends ------------

				if($shorten_array==1)
					{
						
						#print "@current_array\n";	
						my $case_found=0;	
								
						($case_found)=&shorten_array($range,$minimum_no_of_repeats,$sa_allowed_percent_similarity,$user_side_to_shorten,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
						
						

						
						if($case_found==1)
							{
								for(my $m=3;$m>=0;$m--)
									{
										unshift(@modified_array,$current_array[$m]);
										#shift(@modified_array);
									}
								push(@modified_array,$current_array[$#current_array]);	
								
								

								
								@current_array=@modified_array;
								$atleast_one_operation_performed++;
							}
						
						
						undef @modified_array;
					}
				
				
				if($trim_repeats==1 and length($model_repeat)>25)######### $trimming_cutoff set to high ( 20 ) ; $minimum_repeat_length set to 24
					{
						
						
						my @arr_user_sides_to_trim=split(',',$user_side_to_trim);
						
						foreach my $side_and_bases(@arr_user_sides_to_trim)
							{
								
								my $case_found=0;				
								my $tmp_minimum_repeat_length=24;
								my $tmp_trimming_cutoff=20;
								($model_repeat,$case_found)=&trim_repeat_ends_and_improve_score($range,$tmp_minimum_repeat_length,$side_and_bases,$tmp_trimming_cutoff,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
								
								
								if($case_found==1)
									{
										$atleast_one_operation_performed++;
									}
								
								
								if($case_found==1)
									{	
										
										if($side_and_bases=~/NA/)
											{
												undef @modified_array;							
												for(my $m=3;$m>=0;$m--)
													{
														unshift(@modified_array,$current_array[$m]);
													}
												
												my $case_found1=0;
												($case_found1)=&fix_arrays_with_gaps($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
												
												
												push(@modified_array,$current_array[$#current_array]);
																				
												
												if($case_found1==1)	{@current_array=@modified_array;}								
															
												undef @modified_array;	
											}
										#print "\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n";	
										#for(my $k1=0;$k1<=$#original_array;$k1++)
										#	{
										#		print "$original_array[$k1]\n";
										#	} 			
										#print "\n";
										
										#print "Repeat ends trimming operation for $accession:\n\n";
										
										#foreach my $l(@current_array)
										#	{
										#		print "$l\n";
										#	}
										#print "\$model_repeat=$model_repeat\n\n";
										
										#print "\nLF: $left_flank_seq\nRF: $right_flank_seq\n";	
										
										#--- check if the insertion is already handled and removed or not
										#if($model_repeat!~/-/)
										#	{
										#		#$insertion_found_in_repeat=0;
										#	}
								}
									
										
								undef @modified_array;
							}
						
					}
			
			
				if($increase_repeat_length==1 and length($model_repeat)<24)
					{
						#system("echo '$user_side_to_increase_length' >log.txt");
						my @arr_user_side_to_increase_length=split(',',$user_side_to_increase_length);
						
						foreach my $side_and_bases(@arr_user_side_to_increase_length)
							{			
						
								if($side_and_bases!~/NA/)
									{
										my($u_side_to_increase_length,$u_no_of_bases_to_increase_length)=split('-',$side_and_bases);
										#---------------------------------------------------------------
										
										my $case_found=0;	
										undef @modified_array;
										($model_repeat,$case_found)=&increase_all_repeat_lengths($range,$u_side_to_increase_length,$u_no_of_bases_to_increase_length,$accession,$model_repeat,\@current_array,\@modified_array);
										
										
										if($case_found==1)
											{
												#------------ complete the @modified array ---
												for(my $m=3;$m>=0;$m--)
													{
														unshift(@modified_array,$current_array[$m]);
														#shift(@modified_array);
													}
												push(@modified_array,$current_array[$#current_array]);
													
												@current_array=@modified_array;
												
												#undef @modified_array;
												$atleast_one_operation_performed++;
											}
											
										undef @modified_array;	
									}	
							
								else{							#----- auto check with $repeat_extension_identity
								
										#system("echo '$repeat_extension_identity' >log.txt");
										
										my $case_found1=0;	
										my $sides_to_increase_length="NA-0,";
										
										undef @modified_array;
										
										($sides_to_increase_length,$case_found1)=&auto_detect_repeat_sides_and_bases_to_increase($range,$repeat_extension_identity,$accession,$model_repeat,\@current_array,\@modified_array); # special cases: $repeat_extension_identity set to 50
										
										#print "$sides_to_increase_length,$case_found1\n";
										
										if($case_found1==1)
											{
												
												my @arr_sides_to_increase_length=split(',',$sides_to_increase_length);									
												
												foreach my $side_and_bases(@arr_sides_to_increase_length)
													{			
														#system("echo '$side_and_bases' >log.txt");
														
														if($side_and_bases!~/NA/)
															{
																my($side_to_increase_length,$no_of_bases_to_increase_length)=split('-',$side_and_bases);
																#---------------------------------------------------------------
																
																my $case_found=0;	
																
																($model_repeat,$case_found)=&increase_all_repeat_lengths($range,$side_to_increase_length,$no_of_bases_to_increase_length,$accession,$model_repeat,\@current_array,\@modified_array);
																
																
																if($case_found==1)
																	{
																		#------------ complete the @modified array ---
																		for(my $m=3;$m>=0;$m--)
																			{
																				unshift(@modified_array,$current_array[$m]);
																				#shift(@modified_array);
																			}
																		push(@modified_array,$current_array[$#current_array]);
																			
																		@current_array=@modified_array;
																		
																		#undef @modified_array;
																		$atleast_one_operation_performed++;
																	}
																	
																undef @modified_array;	
															}
													}		
											}
											
										undef @modified_array;	
									}	
							
							}	
					
						
						#----- check for misaligned bases
						&check_array_for_misaligned_bases($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
					}			
						
				#--- and re-check the quality score again
				($array_quality_score,$score_det,$score_legend)=&calculate_array_quality_score($range,$species,$accession,$all_gene_positions_folder,$all_gene_positions_file,$matching_reference_repeat,$model_repeat,\@current_array);
				$array_quality_score= sprintf("%.2f",$array_quality_score);
				
			}
			
		
		
		
		
		#################### ask user if they want to see extra degenerated spacer (if present) using 75% identity
		
		if($extend_array==1 and length($model_repeat)>=$word_length)
			{
				
				#@modified_array=@current_array;
				#for(my $m=3;$m>=0;$m--)
				#	{
				#		#unshift(@modified_array,$current_array[$m]);
				#		shift(@modified_array);
				#	}
				#pop(@modified_array);
				
				#print "@current_array\n";	
				#print "Before:\n";
				#foreach my $l(@current_array)
				#	{
				#		print "$l\n";
				#	}
				#print "\t$model_repeat\n";	
				
				
				my $case_found=0;	
						
				($case_found)=&extend_array($range,$max_gap_between_crisprs,$ea_dynamic_search,$ea_allowed_percent_similarity,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
				
				#print "\n******** \$case_found=$case_found *********\n";
				
				if($case_found==1)
					{
						$atleast_one_operation_performed++;
					}
				
				if($case_found==1)
					{
						for(my $m=3;$m>=0;$m--)
							{
								unshift(@modified_array,$current_array[$m]);
								#shift(@modified_array);
							}
						push(@modified_array,$current_array[$#current_array]);	
								
						@current_array=@modified_array;
						undef @modified_array;
						
						#print "After:\n";
						#foreach my $l(@current_array)
						#	{
						#		print "$l\n";
						#	}
						#print "\t$model_repeat\n";	
						
						#------- fix the gaps ---------------------------------
						for(my $m=3;$m>=0;$m--)
							{
								unshift(@modified_array,$current_array[$m]);
							}
						
						my $case_found1=0;
						($case_found1)=&fix_arrays_with_gaps($range,$accession,$model_repeat,$avg_spacer_length,\@current_array,\@modified_array);
						
					
						if($case_found1==1)
							{
								push(@modified_array,$current_array[$#current_array]);
								@current_array=@modified_array;	
								$atleast_one_operation_performed++;
							}
									
						undef @modified_array;	
						
					}
				
				
				undef @modified_array;
			}

					
		#----- skip is length or model_repeat too short, or q_score is too low
		if(($#current_array+1-5) < 3 and (length($model_repeat)< 17 ) or (($#current_array+1-5)<3 and $array_quality_score< -3)  or (length($model_repeat)< 17 and $array_quality_score< -3 ))
			{
				#next;
				$pm->finish;
			}
				
		if($quiet!=1 and $array_quality_score>=0)
			{		
				print "\t  Array_quality_score= $array_quality_score\n";
			}
			

			
		#----------------- Very Important: save a copy of the @current_array : will be used for gff creation ----------
		my @backup_current_array=@current_array;
		
			
		#----------- reverse the array ---------------------------------------------------------
		
		if($user_reverse_the_array==1){$array_direction = "R";}
		if($array_direction =~ /R/)
			{
				my $case_found=0;	
				($model_repeat,$case_found)=&revers_an_array($range,$accession,$model_repeat,\@current_array,\@modified_array);
				
				
				if($case_found==1)
					{
						if($potential_alternate_repeat ne "NA")
							{
								$potential_alternate_repeat=reverse $potential_alternate_repeat; $potential_alternate_repeat=~tr/ACGT/TGCA/;
							}
							
						
						#------------ complete the @modified array ---
						for(my $m=3;$m>=0;$m--)
							{
								unshift(@modified_array,$current_array[$m]);
								#shift(@modified_array);
							}
						push(@modified_array,$current_array[$#current_array]);
							
						@current_array=@modified_array;
						
						#undef @modified_array;	
						$atleast_one_operation_performed++;
					}
				undef @modified_array;	
			}
		



		
		



		#---------------------------------------------------------------------------------------------------------------------------------------------------------------
				
		
		
		#---- finally finalize the array with fixing Deletion texts --------------------
		my $finalize_array=1;
		my $questionable_array=0;
		my $array_start_position=0;
		my $array_stop_position=0;
		
		if($finalize_array==1)#if($gap_found_in_repeat==1 or $insertion_found_in_repeat==1)
			{
				#print "\nA1: $left_flank_seq\nA2: $right_flank_seq\n\n";
				
				my $case_found=0;

				#undef @backup_current_array;
				#@backup_current_array=@current_array;
					
				($model_repeat,$potential_alternate_repeat,$questionable_array,$array_start_position,$array_stop_position,$case_found)=&finalize_the_array($range,$array_quality_score,$score_det,$score_legend,$all_gene_positions_folder,$all_gene_positions_file,$array_direction,$matching_reference_repeat,$repeat_family,$array_direction_MEMO,$potential_alternate_repeat,$accession,$model_repeat,$avg_spacer_length,$left_flank_length,$right_flank_length,\@current_array,\@modified_array);
				
				
				if($case_found==1)
					{
						$atleast_one_operation_performed++;
					}								
			}
		
		#print "All done\n";
		
		
		#------------------------------ now store the array if a copy of it is not already stored -------------------------------------------------
		
		my $skip_c2=0;
		if($skip_c2==1)
			{
				my $array_already_exist=0;
				
				my $a_s_p;my $a_st_p;
				if($array_start_position<$array_stop_position)
					{	 
						$a_s_p=$array_start_position;
						$a_st_p=$array_stop_position;
					}
				else{
						$a_s_p=$array_stop_position;
						$a_st_p=$array_start_position;
					}		
				
				
				
				#print "\n\n";
				#foreach my $q_score(sort{$b<=>$a} keys %hash_of_arrays_per_accession)
				#  {
					foreach my $asp(sort{$a<=>$b} keys %hash_of_arrays_per_accession)
					 {
						#print "\$asp=$asp\n";
						foreach my $astp(sort{$a<=>$b} keys %{$hash_of_arrays_per_accession{$asp}})
							{	
								#print "\$astp=$astp\n";
								
								my $lb; my $ub;
								if($asp<$astp){$lb=$asp;$ub=$astp;}
								else{$lb=$astp;$ub=$asp;}
								
								my $array_center=int($lb+($ub-$lb)/2);
								my $existing_middle_point=int($lb+($ub-$lb)/2);
								my $current_middle_point=int($a_s_p+($a_st_p-$a_s_p)/2);
								
								#print "START: $asp\tSTOP: $astp\t LB:$lb\tUB:$ub \t\$array_start_position=$array_start_position\t\$array_stop_position=$array_stop_position\n";
								
								#---- check if either $array_start_position or $array_stop_position matches with existing record or not --------------------------------------
								
								if($a_s_p==$lb or $a_st_p==$ub or ($current_middle_point>$lb and $current_middle_point<$ub) or ($existing_middle_point>$a_s_p and $existing_middle_point<$a_st_p))
									{
										#----- check which one is longest, and delete the other one is shorter than the current one ---------------------------------------
										
										
										if((abs($ub-$lb) < abs($a_st_p-$a_s_p)) and ($existing_middle_point>$a_s_p and $existing_middle_point<$a_st_p))
											{
												#print "\tExisting array shorter than current one\n\n";
												delete $hash_of_arrays_per_accession{$asp}{$astp};
											}
										elsif($current_middle_point>$lb and $current_middle_point<$ub)
											{
												$array_already_exist=1;
												last;
											}	
									}
								#elsif($array_center>$a_s_p and $array_center<$a_st_p)
								#	{
								#		$array_already_exist=1;
								#		last;
								#	}								
							}
						if($array_already_exist==1){last;}	
					}
				 # }
			}
		
		
		my ($array_already_exist,$existing_array_q_score)=&check_array_existance($array_quality_score,$array_start_position,$array_stop_position,$tmp_output_file);
		if($array_already_exist==1)
			{
				if($quiet!=1 and $array_quality_score>=0)
					{
						#print "\t\tAlready part of previously processed Array. Skipping...\n";
					}	
				#next;
				$pm->finish; 
			}
		
			
		if($array_already_exist==0 and $array_quality_score>=$array_quality_score_cutoff and length($model_repeat)>=$repeat_length_cutoff and ($#backup_current_array+1-5) >= $minimum_no_of_repeats)
			{	
				#--------- append the positions of spacer and repeats in the GFF file ------------------------------------------------------------------------------------------
				my($array_start_pos1,$array_stop_pos1,$avg_spacer_len1)=&find_array_start_stop_position('F',\@backup_current_array);			
				&create_gff_file($gff_file,$array_direction,$crispr_index,$array_start_pos1,$array_stop_pos1,$accession,$model_repeat,\@backup_current_array,\@modified_array,\%hash_id_lookup_table);
				
				
				my $current_array_content=join("&&&",@current_array);
				
				open(APP,">>$tmp_dir\/$tmp_output_file") or print "$!";
				flock(APP,2);
				print APP "$accession%%%$array_start_position%%%$array_stop_position%%%$array_quality_score%%%$current_array_content\n";
				close(APP);

							
				$crispr_index++;
			}			
		elsif($array_already_exist==0 and $array_quality_score>=0 and length($model_repeat)>=$repeat_length_cutoff)
			{
				#----- dump the poor arrays in junkfile
				
				open(FP,">>$filtered_out_crisprs");
				flock(FP,2);
				my $fp_array=join("\n",@current_array);
				if(defined $hash_id_lookup_table{$accession})
				{
					$fp_array=~s/($accession)\|?/$hash_id_lookup_table{$accession}/g;
				}	
				print FP "$fp_array\n\n";
				close(FP);
			}

				
		#---------------------------------------------------------------------------------	

	}
		
		undef @current_array;	
		
		$pm->finish;
	
	}


	if($continue_from_last_process ==1)
		{
			#---- add the current accession to continue file -------------------
			open(ADD,">>$tmp_dir\/$continue_from_last_process_file");
			flock(ADD,2);
			print ADD "$accession\n";
			close(ADD);
			#-------------------------------------------------------------------
		}
	
	
	$pm->wait_all_children;



 
	#-------------------- now load the arrays in %hash_of_arrays_per_accession -----------------------------------------------------
	open(RD4,"$tmp_dir\/$tmp_output_file") or print "$!";
	flock(RD4,2);
	my @arr_rd4=<RD4>;
	close(RD4);
	foreach my $line(@arr_rd4)
		{
			my @arr_line=split('%%%',$line);
			my $accession=$arr_line[0];
			my $start=$arr_line[1];
			my $stop=$arr_line[2];
			my $q_score=$arr_line[3];
			my $existing_array=$arr_line[4];
			
			if(defined $start and defined $stop and defined $existing_array)
				{			
					$hash_of_arrays_per_accession{$start}{$stop}=$existing_array;						
				}	
		}
	


	#-------------------- now print the array(s) for the current accession accession -----------------------------------------------
	#print "\n\nCRISPRDetect output:\n\n";	
	
	
	open(WR,">>$output_file");
	flock(WR,2);
	my $array_index=1;
	
	#foreach my $q_score(sort{$b<=>$a} keys %hash_of_arrays_per_accession)
	#	{			
			foreach my $start_pos(sort{$a<=>$b} keys %hash_of_arrays_per_accession)
				{
					foreach my $stop_pos(sort{$a<=>$b} keys %{$hash_of_arrays_per_accession{$start_pos}})
						{
							#print "$hash_of_arrays_per_accession{$start_pos}{$stop_pos}\n\n";
							print WR "Array $array_index $start_pos-$stop_pos \t\t**** Predicted by CRISPRDetect $version *** \n";
							my $cd_array=$hash_of_arrays_per_accession{$start_pos}{$stop_pos};
							if(defined $hash_id_lookup_table{$accession})
							{
								$cd_array=~s/($accession)\|?/$hash_id_lookup_table{$accession}/g;
							}	
							
							$cd_array=~s/&&&/\n/g;
							#print "$1\n";
							print WR "$cd_array\n\n";
							$array_index++
						}
				}
	#	}	
	close(WR);
	
	unlink("$tmp_dir\/$tmp_output_file");		
	
	
	
	#------- delete files not needed anymore ----
	if(-e "$tmp_dir\/$accession\.fna")
		{
			if($keep_input_sequences==0)
			{
				unlink("$tmp_dir\/$accession\.fna");
			}	
		}		
	#unlink("$tmp_dir\/$accession\.fna");	
	
	#--------------------------------------------
	
	#$pm->finish;
	
	#sleep(1);
 }
 
#$pm->wait_all_children;

#----- clean $tmp_dir\/ folder --------------- 
#unlink<$tmp_dir\/*\.mfe>;
#unlink<$tmp_dir\/mr_*\.txt>;

#print "\n\nall_gene_positions_file = $tmp_dir\/$all_gene_positions_file\n\n";
#unlink("$tmp_dir\/$all_gene_positions_file");

if($blast_db_file_of_known_repeats=~/^$tmp_dir/)
	{
		#print "\n\nRemoving all files associated to $blast_db_file_of_known_repeats\n\n";
		system("rm $blast_db_file_of_known_repeats\*");
		
	}

foreach my $seq_file(@arr_sequence_files)
	{
		#print "\$seq_file=$seq_file\n";
		if($keep_input_sequences==0)
			{
				unlink("$tmp_dir\/$seq_file\.fna");
			}			
	}
#print "\nAll done. :)\n\n";








##################################################################### subs #############################################################

sub check_executables_required_for_CRISPRDetect()
	{

		#---------------- check CLUSTALW installation -----------------------------------------------------------
		my $clustalw_path=`which clustalw >&1 2>&1`; chomp $clustalw_path; $clustalw_path=~s/\r//g;

		my $clustalw_found="Found in $clustalw_path";

		if(not defined $clustalw_path or $clustalw_path eq "" or $clustalw_path=~/:/)
		{
			if(-e 'clustalw')
				{
					my @clustal_help=`clustalw -help >&1`;
					my $clustalw_identified=0;
					foreach (@clustal_help)
						{
							if($_=~/CLUSTAL (\d)\./)
								{
									#print "$_\n";	
									if($1>=2)
										{						
											$clustalw_identified++;
										}
									else{								  
											print qq~\nNOTE: Although clustalw was found in current directory, the current version is lower than Version 2. Please download and install clustalw [ http://www.clustal.org/clustal2/ ].\nOnce you download and install clustalw2, rename it to clustalw.~;
										}		
									last;
								}
						}
						
					if($clustalw_identified>0)
						{
							$clustalw_found="Found in current directory";
						}
					else{
							$clustalw_found="Not found";
						}		
				}
			else{
					print qq~\nNOTE:\nclustalw was not found in system path. Please download and install clustalw [ http://www.clustal.org/clustal2/ ]. Version 2 is required. Once you download and install clustalw2, rename it to clustalw.~;
					$clustalw_found="Not found.";
				}
		}

		#---------------- check RNAfold installation ------------------------------------------------------------

		my $RNAfold=`which RNAfold >&1 2>&1`; chomp $RNAfold; $RNAfold=~s/\r//g;

		my $RNAfold_found="Found in $RNAfold";

		if(not defined $RNAfold or $RNAfold eq "" or $RNAfold=~/:/)
		{
		print qq~   
		   NOTE:
			RNAfold was not found in system path. Please download and install RNAfold (version 2.0.7 or higher)
			from http://www.tbi.univie.ac.at/\~ronny/RNA/index.html#download
			Latest version [V2.1.3]: [ http://www.tbi.univie.ac.at/\~ronny/RNA/packages/source/ViennaRNA-2.1.3.tar.gz ].
			
			The related method 'Analysis of RNA secondary structure' will not be available without RNAfold.
			
		~;

			#$check_secondary_structure_of_repeat=0; #-- skipping this method as RNAfold is not available
			$RNAfold_found="Not found.";
		}
		else{	
				#--- check version >= 2.0.7 ----------------------------------------
				my $RNAfold_version=`$RNAfold -V >&1 2>&1`;chomp $RNAfold_version; $RNAfold_version=~s/\r//g;
				
				$RNAfold_version=~s/RNAfold //; $RNAfold_version=~s/\.//g;
				
				if(not $RNAfold_version>=207)
					{
						#$check_secondary_structure_of_repeat=0; #-- skipping this method as RNAfold is not available
						$RNAfold_found="Older version found.";
						
						print "\n    Note: It seems that you have an older version of RNAfold, please update your version from http://www.tbi.univie.ac.at/~ronny/RNA/index.html#download\n\n";
					}	
		}


		#------------ Check seqret installation----------------------------------------------------------------------------

		my $seqret_path=`which seqret >&1 2>&1`; chomp $seqret_path; $seqret_path=~s/\r//g;

		my $seqret_found="Found in $seqret_path";

		if(not defined $seqret_path or $seqret_path eq "" or $seqret_path=~/:/)
		{
			if(-e 'seqret')
				{
					my @seqret_help=`seqret -h >&1`;
					my $seqret_identified=0;
					foreach (@seqret_help)
						{
							if($_=~/Version:/)
								{
									$seqret_identified++;
									last;
								}
						}
					if($seqret_identified>0)
						{	
							$seqret_found="Found in current directory.";
						}
					else{
							$seqret_found="Can not execute seqret present in current directory.";
						}		
				}
			else{
			print qq~
			NOTE: 
			  seqret was not found in system path. Please download and install seqret [ http://emboss.sourceforge.net/download/ ]
			  version EMBOSS:6.3.1 tools or higher required.	
			
			  'seqret' is required for extracting genomic sequence from GBK files.
			~;
			$seqret_found="Not found.";
			}
		}
		#--------------------------------------------------------------------------------------------------------------------

		#------------ Check water installation----------------------------------------------------------------------------

		my $water_path=`which water >&1 2>&1`; chomp $water_path; $water_path=~s/\r//g;

		my $water_found="Found in $water_path";

		if(not defined $water_path or $water_path eq "" or $water_path=~/:/)
		{
			if(-e 'water')
				{
					my @water_help=`water -h >&1`;
					my $water_identified=0;
					foreach (@water_help)
						{
							if($_=~/Version:/)
								{
									$water_identified++;
									last;
								}
						}
					if($water_identified>0)
						{	
							$water_found="Found in current directory.";
						}
					else{
							$water_found="Can not execute seqret present in current directory.";
						}
				}
			else{	
			print qq~			
			NOTE: 
			   water was not found in system path. Please download and install water [ http://emboss.sourceforge.net/download/ ]
			   version EMBOSS:6.3.1 tools or higher required.	
			
			   'water' is required for comparing reference known repeats and CRISPR array repeats.
			~;
			
			$water_found="Not found.";
			}
		}
		#--------------------------------------------------------------------------------------------------------------------


		if($RNAfold_found=~/Not/ and $water_found=~/Not/)
			{
				
		print qq~		
		   NOTE:
			You need to install the Emboss tools (version 6.3.1 or higher) and Vienna RNA package version 2.0.7 or higher. The other_executables folder
			has only the executables of RNAfold, water and seqret. Check each of them by the following commands from the current directory:
			
			$cd_path/other_executables/RNAfold -h
			
			$cd_path/other_executables/water -h
			
			$cd_path/other_executables/seqret -h
			
			If you see all the programs are running successfully, then copying them to a location like '/usr/local/bin' [require root access] will solve the issue.
			Alternately export the 'other_executables' folder to your path will also work. To know how to export a directory refer this http://www.troubleshooters.com/linux/prepostpath.htm
		~;		
			}

		
		#------------ Check cd-hit-est installation----------------------------------------------------------------------------

		my $cd_hit_path=`which cd-hit-est >&1 2>&1`; chomp $cd_hit_path; $cd_hit_path=~s/\r//g;

		my $cd_hit_found="Found in $cd_hit_path";

		if(not defined $cd_hit_path or $cd_hit_path eq "" or $cd_hit_path=~/:/)
			{

			print qq~
			NOTE:
			   cd-hit-est was not found in system path. Please download and install cd-hit-est [ http://weizhongli-lab.org/cd-hit/download.php ]
			~;
			$cd_hit_found="Not found.";
			}
		
		#------------ Check blastn installation----------------------------------------------------------------------------

		my $blastn_path=`which blastn >&1 2>&1`; chomp $blastn_path; $blastn_path=~s/\r//g;

		my $blastn_found="Found in $blastn_path";

		if(not defined $blastn_path or $blastn_path eq "" or $blastn_path=~/:/)
			{

			print qq~
			NOTE:
			   blastn: 2.2.30+ was not found in system path. Please download and install blastn [ ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ ]
			~;
			$blastn_found="Not found.";
			}
		#--------------------------------------------------------------------------------------------------------------------		
		
		return ($clustalw_found,$water_found,$seqret_found,$RNAfold_found,$cd_hit_found,$blastn_found);
}


sub show_help()
	{	
			
		my($clustalw_found,$water_found,$seqret_found,$RNAfold_found,$cd_hit_found,$blastn_found)=@_;
		

print qq~
CRISPRDetect Version:2.2 help:

Syntax:	
	perl CRISPRDetect.pl -g NZ_CP006019.gbk -o NZ_CP006019_CRISPRDetect > NC_003106_CRISPRDetect.log
	
	     Runs with defaults on a complete gbk file that has cas1 or cas2 annotated
	
	perl CRISPRDetect.pl -f test_multifasta.fa -o test_CRISPRDetect -check_direction 0 -array_quality_score_cutoff  3 -T 0 > test.log
        
        Runs with a lower score cutoff on a fasta file, cutoff 3 rather than 4, as cas1 and cas2 are not annotated and would score +1. Appropriate for contigs/fasta. -T 0, use all processors rather than the default of 4. Does not check direction (not recommended) 
	
~;

print color('bold blue');
print qq~

CRISPRDetect dependencies:
	clustalw 	$clustalw_found.	
	water 		$water_found 		[Comes with EMBOSS:6.3.1 tools]
	seqret 		$seqret_found 		[Comes with EMBOSS:6.3.1 tools]
	RNAfold 	$RNAfold_found		[Comes with Vienna RNA package]
	cd-hit-est 	$cd_hit_found 	[Comes with cdhit package] 	
	blastn 		$blastn_found 		[Comes with ncbi-blast+ package] 
~;
print color('reset');
printf qq~




Compulsory parameters:
---------------------
 
	-f/-g		FASTA/Genbank	Specify a FASTA formatted [e.g. -f test.fa] or Genbank formatted file [e.g. -g NZ_CP006019.gbk] containing the sequence. 		
	
	Note: the default cutoff of 4 is appropriate for genbank files that have cas1 or cas2 annoated, 3 is more appropriate for fa.
	
	-o		TEXT		Specify a text file that will contain the output [e.g. -o NC_003106_CRISPRDetect] 			
						Note: CRISPRDetect will provide two additional output files - one containing the filtered out arrays (e.g. NC_003106_CRISPRDetect.fp) 
						and a gff annotation file (e.g. NC_003106_CRISPRDetect.gff)
						
Basic options:
-------------						
	-h/-help	HELP		shows this help text
	-q/-quiet	0 or 1		Switch off/on step by step progress reporting [default 0]	
	-T		Threads		Specify number of parallel processes CRISPRDetect should use  [default 4; specify '-T 0' to use all processors]		
	-tmp_dir	tmp/		This is the default directory where temporary files generated by CRISPRDetect and its dependencies will be stored.	


Parameters for putative CRISPR identification [optional]:
--------------------------------------------------------	
	-word_length			11	This is the default word length CRISPRDetect uses to find the putative CRISPRs. Any positive integer >=6 can be used.
	-minimum_word_repeatation	3	By default CRISPRDetect uses 3 repeating identical words to find putative CRISPRs. To find CRISPRs with 2 repeats, use -minimum_word_repeatation 2	
	-max_gap_between_crisprs	125	By default the maximum gap is set tp 125 nucleotides between the repeating identical seed words.
	-repeat_length_cutoff		17	After the intial processing, putative CRISPRs with repeat lengths less than this value will be rejected.


Filtering parameters [optional]:
-------------------------------	
	-minimum_repeat_length		23	Minimum length of repeats 
	-minimum_no_of_repeats		3	Predicted CRISPRs with number of repeats less than this value will be excluded. To include CRISPRs with only 2 repeats, use -minimum_no_of_repeats 2
	-array_quality_score_cutoff	4	Predicted CRISPRs with score less than this value will be excluded from the output file. 
						The CRISPRs with score >=0 and less than the specied value will be moved to the output.fp file [output refers to user given output filename]. Cutoff of 3 is more appropriate for fasta files.
						

Additional parameters [optional]:
--------------------------------
	-left_flank_length		500	This is the default length of the 5' (upstream) region of the CRISPRs.
	-right_flank_length		500	This is the default length of the 3' (downstream) region of the CRISPRs.		
	


Advanced options [optional]:
---------------------------	
	 
	To test different methods as specified in the literature, open the CRISPRDetect.pl program with any text editor [e.g. gedit in RHEL/Fedora/CentOS, or vi in any *nix OS, or 
	notepad in Windows OS] and change the parameters in the top most section of the script. To toggle individual methods, locate the '\$check_' prefix and change the value to 1 
	(i.e. the method will be applied) or 0 (i.e. the method will not be applied). 
	
	Examples:
		
		Direction specific options:
		--------------------------
			\$check_direction=0;			[ Default is 1, making it 0 will turn the method off.] 
				
		To change the parameter(s) of a particular method (such as check_array_degeneracy) change the nested variables under that particular method.
		
			\$check_array_degeneracy=1;	 
				\$array_degeneracy_score=0.41; 		[ Default: PPV (0.91) - 0.50 ]
				\$permitted_mutation_per_array=0; 	[ Default 0 ]
	
	Changing to '\$permitted_mutation_per_array=2;' will instruct the program to allow maximum 2 bases as permitted mutations per CRISPR array.




NOTE:   
	Please make sure that the 'clustalw', 'RNAfold', 'water' , 'seqret', 'cd-hit-est' and 'blastn' are in the system path and have execution permission for the current user. The 'tmp'
	folder in the current directory should have read and write permissions. An easy way to do that is by issuing the command 'chmod -R 755 . && chmod 777 tmp' from the current directory.

	CRISPRDetect.pl should run under any unix based operating system that has a working 'perl' executable [comes with default installations under all *nix based operating systems]. 
	However, Mac OS users may require to  compile and install RNAfold (comes with vienna RNA package).
           
	For version updates and bug fixes refer to http://bioanalysis.otago.ac.nz/CRISPRDetect  

 	
~;

return 1;
}

sub get_path()
	{
		my $dir=shift(@_);
		my @arr_p1=split('\/',$cd_path);
		pop(@arr_p1);
		$dir=join("\/",@arr_p1);
		
		return $dir;
	}







############################################################### end of subs #################################################################	
