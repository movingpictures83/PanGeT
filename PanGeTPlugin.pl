#!/usr/bin/perl
use warnings;
use strict;

my $inputfile;

sub input {
   $inputfile = $_[0];
}

sub run {
}

sub output {

use Cwd qw();
my $path = Cwd::cwd();
print "$path\n";

open (PARA, $inputfile) or die "FAILURE: $inputfile: $!";

my $pathto; 
my $out_path; 
my $hcutoff; 
my $e_value;

my $no_of_threads;
my $species;


my @genome_list;
my @parameter=<PARA>; 

my $t_gene='-';

my $flag = 0;

my %hash_name; 

my $id; 

my $hvalue;


my $h;
my %hvalue;
foreach my $ip(@parameter)
{

	if($ip=~/Input_directory\s+(.+)/)
	{
	$pathto=$1;


	($species) = $ip=~m/.+\/(\w+)/g;
	}
	if($ip=~/Output_directory\s+(.+)/)
	{
	 $out_path=$1;
	}
	if($ip=~/Hvalue\s+(.+)/)
	{
	$hcutoff=$1;
	}
	if($ip=~/Evalue\s+(.+)/)
	{
	 $e_value=$1;
	}
	if($ip=~/Threads\s+(.+)/)
	{
	 $no_of_threads=$1;
	}
}
print "PanGET\t\t\t\n";
print "Input path= $pathto\tOutput path=$out_path\tH value=$hcutoff\tEvalue=$e_value\tThread=$no_of_threads";
                             
chdir("$pathto");              #entering the selected directory
print"\n\n";


my $files=`ls *.fna`||die ("\"please check the path of input directory!!!!!!!!!\"");     #list all the .fna file based on that creating one genome compared file

$files=~s/\.fna//g;


open(GENOME,">genome_compared")||die("couldnt open file");
print GENOME "$files";

my $count_unique=0;          #intialization
my $count_core=0;
 
open(GENOME,"genome_compared")||die("couldnt open file");
@genome_list=<GENOME>;                             #list of the strains in the genome_list array
chomp(@genome_list);
my @whole=@genome_list;

#print "enter the Number of the strain you want keep as reference\n\n";

#for my $th(0..$#genome_list)
#{
#	print "$th\.$genome_list[$th]\n";                            # printing the folder name with numbers

#}

my $va=0;#<STDIN>;
chomp($va);

my $first=$genome_list[$va];    
chomp($first);
                      
print "you have taken: $first as reference strain \n\n\n";                              #getting the reference strain
close GENOME;












`mkdir -p $out_path/process_$species`;                  #creating folders in output_files
`mkdir -p $out_path/output_$species`;
`mkdir -p $out_path/process_$species/common`;


my @name1;
my $count_genome;
my $number;
my @code_length;
open(NUMBER,">$out_path/output_$species/no_of_genes.txt")||die("couldnt open file");  

foreach my $element(@genome_list)
{

		$element=~s/\s//g;                      #count the no the strains leaving the empty line given the genome compared

		if ($element=~/^$/)
		{
		next;
		}
		if ($element!~/^$/)
		{
		$number++;
		$count_genome++;
		}




chomp($element);
open(FILE,"$element.fna")||die("couldnt open file $element.fna");       #opening the each strain .fna in a loop
my $first_line = <FILE>;

my $genome_seq=do{local $/; <FILE>};
close FILE;

push( @name1,$first_line);               

# delete the header in the fasta file

$genome_seq=~s/\n//g;


                                   ###########    CDS FILE_CREATION    #############
open(EDIT,">$out_path/process_$species/common/$element.edit.ptt")||die("couldnt open file");    #rewriting the .ptt file as edit.ptt with only needed information
open(CDS,">$out_path/process_$species/common/$element.cds")||die("couldnt open file");	#seqfile 
open(PPT,"<$element.ptt")||die("couldnt open file");


#print CDS "$name1\n";	


my $count_line;
my $maxcode_len;
my $len_code;
my $portion;
my $code;

my $PID;
my $por;
my $max;
my $product;	

my $geneid;


my $COG;
while(<PPT>)
{
		
		if($_=~/^(\d+)\.\.(\d+)\s*(\W)\s\d*\s*(\d*)\s*(\S*)\s(\S*)\s*\W\s*(\S*)\s*(.*)/)

		{
			$count_line++;
			my $start_CDS=$1;
		my $end_CDS=$2;
		my $strand_symbol=$3;
		$PID=$4;
		$geneid=$5;
		$code=$6;
		$COG=$7;
		$product=$8;
        	my $len_code=length($code);
        if($count_line==1)
		{		
        $maxcode_len=$len_code;
		}
        if($len_code>$maxcode_len)
		{
		 $maxcode_len=$len_code;
		}
		print EDIT "$start_CDS..$end_CDS\t$strand_symbol\t$code\n\n";
			
                       	
		  
                if($strand_symbol eq '+')
                {
                my $portion2=substr($genome_seq,$start_CDS-1,($end_CDS-$start_CDS+1));
                print CDS">$start_CDS..$end_CDS:$code\t\t$strand_symbol\t\t$geneid\t\t<a href=http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=kegg&dbkey=kegg&keywords=$code>$code</a>\t\t$COG\t\t$product\n$portion2\n";
                }
                else
                {
                my $por1=substr($genome_seq,$start_CDS-1,($end_CDS-$start_CDS+1));
                my $portion3=reverse($por1);
                $portion3=~tr/ATGC/TACG/;
                print CDS">$start_CDS..$end_CDS:$code\t\t$strand_symbol\t\t$geneid\t\t<a href=http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=kegg&dbkey=kegg&keywords=$code>$code</a>\t\t$COG\t\t$product\n$portion3\n";
                }
            


		}
         

	}
	print NUMBER "$element=$count_line\n";
	

    push(@code_length,$maxcode_len);
	
}



@code_length=sort {$a <=> $b} @code_length;
my $max=$code_length[-1];
 


print "CDS FILE CREATED\n";
close EDIT;
close CDS;
close PPT;
close NUMBER;



while ($first ne $genome_list[0])
{
	my $b=shift(@genome_list);                         #shuffling the genome list

	push(@genome_list,$b);
}



############################# whole running block ########################################################


my @hvalue;
my @geneid;
my @seq;
my $flag_val;
my %id;
my %h_value;

my $no;
my $name;
for my $y(0..$#genome_list)
{
		$no++;
		if($no>1 && $no<$#genome_list+2)
		{
		$b=shift(@genome_list);                         #shuffling the genome list

		push(@genome_list,$b);
		}
		$name=$genome_list[0];
		chomp($name);

		`mkdir -p $out_path/process_$species/$name`;



		if($y==0)
		{
		open(OUT,">$out_path/output_$species/name")||die("couldntopenfile");
		print OUT "@name1\n";
		close OUT;
		}

		open(OPE,"$out_path/process_$species/common/$name.edit.ptt")||die("couldnt open file");
		open(WRITE,">$out_path/process_$species/$name/list")||die("couldnt open file");	
		while	(<OPE>)
		{
			if($_=~/^\d+\.\.\d+\s+[+-]\s+(.+)/)#taking the Synonym of edit.ptt
			{
			 my $ref=$1;
			
			print WRITE"$ref\n";               #taking a list of all genes and write in list 
			}
		}
		close OPE;
		close WRITE;
	
		print "\n Reference genome =$genome_list[0]\n";

		for (my $i=1;$i<=$#genome_list;$i++)
		{
			chomp($genome_list[0]);
			chomp($genome_list[$i]);
			my $taken_genome=$genome_list[0];
			my $next_genome=$genome_list[$i];
			my $output="$taken_genome.$next_genome.result";
			
			my $hash_name="$taken_genome\_$next_genome";
							            ##### RUNNING BLAST ####


			`makeblastdb -in $next_genome.fna -dbtype 'nucl' `;
			`blastn -query $out_path/process_$species/common/$taken_genome.cds -db $next_genome.fna  -num_descriptions 1 -num_alignments 1 -evalue $e_value -num_threads  $no_of_threads  -out $out_path/process_$species/$name/$output `;



			print "$i blast over\n";

			my $a=$i+1;
			open(BLAST1,">$out_path/process_$species/$name/1\_$a")||die("couldnt open the file");

			open(RESULT,"$out_path/process_$species/$name/$output")||die("couldnt open the file");
			my @array=<RESULT>;
			close RESULT;
			open(GD,"$out_path/process_$species/common/$next_genome.edit.ptt")||die("couldnt open the file");
			my @next=<GD>;
			close GD;
 #### READING BLAST OUTPUT FILE###
				                   #this part varies based on the blast version 
			my @querystart;
			my @subject_start;                    # Array of starting and ending positions in the alignment
			my @queryend;
			my @subject_end;
			my @Identity;
			my @subject_seq;


			my $startmatch;
			my $start_match;

			my $start_query;
			my $end_query;

			my $st;
			my $en;
			my $id1;



			my ($len_query, $length_match,$percen,$score,$hvalue,$target_strand,$query_strand);

				for my $g(0..$#array)
				{

					if($array[$g]=~/BLASTN/)               
					{
				       
					}

					if($array[$g]=~/^Query\=\s+c?\d+\.\.\d+\:(.+)/) # Maching query,in the fifth line NO of hits (counting no of uniques) 
					{
				   	
						$score=0;           #intializing score =0;
						 $flag=0;	
						 $id=$1;
	
						$hvalue='0.0';
						$t_gene='-';
	
						if($array[$g+2]=~/^Length=(\d+)/)  #taking length 
						{
							$len_query=$1;
						}
						if($array[$g+5]=~/\*+\s+No\s+hits\s+found\s+\*+/)
						{
							print BLAST1 "$id\t$hvalue\t-\n"; 
							$len_query='';
							$id='-';
						}

					}




					if($array[$g]=~/\s+Identities\s+\=\s+\d+\/(\d+)\s+\((\d+)\%\)/)
					{
	
					 $length_match=$1;                      #Taking Identity 
					 $percen=$2;
					}

					if($array[$g]=~/\s?Strand\s?\=\s?(\w+)\s?\/\s?(\w+)/) #Strand=Plus/Plus
					{
						  $query_strand=$1;
						 $target_strand=$2;  

						$target_strand=~s/Plus/\+/;
						$target_strand=~s/Minus/\-/;     

					}
	
					if($array[$g]=~/^\s+Score\s+\=\s+/)
					{
						$score++; 
		                             #counting the scores
						@querystart=();
						@subject_start=();                    # Array of starting and ending positions in the alignment
						@queryend=();
						@subject_end=();
						@Identity=();


					       
					}
	

	

	
					if($array[$g]=~/^Query\s+(\d+)\s+[A-Za-z\-]+\s+(\d+)/)
					{
						
						 
						$start_query=$1;
	
						push(@querystart,$1);             
						$end_query=$2;

						push(@queryend,$end_query);
						$flag=$g+4;
				      
					}
	
					if($array[$flag]!~/^Query / && $flag!=0 && $flag==$g)
					{
						if($score==1 )
						{
	
							 $hvalue=(($percen/100)*$length_match)/$len_query;
	
						 	$hvalue=substr($hvalue,0,4);
						



						$startmatch=substr($subject_start[0],0,3);
						$start_match=substr($subject_end[$#subject_end],0,3);

					       
							foreach (@next)
							{	
								if($_=~/$startmatch/||$_=~/$start_match/)
								{
	
									if($_=~/^(\d+)\.\.(\d+)\s+[\+-]\s+(\w+\.?\_?\w?)/)
									{
										$st=$1;
										$en=$2;
										$id1=$3;
									

										 $t_gene='-';
										
										if($target_strand eq '+')
										{
											if($subject_start[0]>=($st-200) && $subject_start[0]<=($st+200))
											{
												$t_gene=$id1;
												last;
											}
				
										}
										else
										{
											if($subject_end[$#subject_end]>=($st-200) && $subject_end[$#subject_end]<=($st+200))
											{
												$t_gene=$id1;
												last;
											}
										}
			

											

									}	
								}
							}

					#	print "$t_gene\n";
			
						if(($t_gene eq '-')&& ($name eq  $first))
						{
						my $start_CDS;
						my $end_CDS;
						my $strand_symbol;
						my $geneid = "";
						my $code = "";
						my $COG = "";
						my $product ="";
						my $por1;
						
						my $portion2;
						my $portion3;
						my $start;
						my $end;
						open(VG,">>$out_path/process_$species/common/$name.edit.ptt")||die("couldnt open file");	

						print VG "$subject_start[0]..$subject_end[$#subject_end]\t$target_strand\tUTR$subject_start[0]\_$subject_end[$#subject_end]\n\n";
						my $subj_seq=join('',@subject_seq);
						
						 $t_gene='UTR_'.$subject_start[0].'_'.$subject_end[$#subject_end];
						
						open(CD,">>$out_path/process_$species/common/$next_genome.cds")||die("couldnt open file $next_genome.cds"); 
						open(FQ,"$next_genome.fna")||die("couldnt open file $next_genome.fna");  
						

						my $first_line = <FQ>;	
						push(@name1,$first_line);
               					my $genome_seq = do{local $/; <FQ> };

						$genome_seq=~s/\n//g;
						$start_CDS=$subject_start[0];
						$end_CDS=$subject_end[$#subject_end];	
						$strand_symbol=$target_strand;
						my $PID=$t_gene;
				
					     
						if($strand_symbol eq '+')
						{
						 $portion2=substr($genome_seq,$start_CDS-1,($end_CDS-$start_CDS+1));
						 print CD ">$start_CDS..$end_CDS:$PID\t\t$strand_symbol\t\t$geneid\t\t<a href=http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=kegg&dbkey=kegg&keywords=$code>$code</a>\t\t$COG\t\t$product\n$portion2\n";
						}
						else
						{
						my $por1=substr($genome_seq,$start_CDS-1,($start_CDS-$end_CDS+1));
						$portion3=reverse($por1);
						$portion3=~tr/ATGC/TACG/;
					
						print CD">$start_CDS..$end_CDS:$PID\t\t$strand_symbol\t\t$geneid\t\t<a href=http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=kegg&dbkey=kegg&keywords=$code>$code</a>\t\t$COG\t\t$product\n$portion3\n";
						}
					   
				
						
 						}



	
					close CD;
						close FQ;
						close VG;
						{no strict "refs";

						if($hvalue>$hcutoff)
						{
							$$hash_name{$id} = $t_gene;
						}

				}
						print BLAST1 "$id\t$hvalue\t$t_gene\n"; 
						@querystart=();
						@subject_start=();                    # Array of starting and ending positions in the alignment
						@queryend=();
						@subject_end=();
						@Identity=();
						
					 	
	
						}
					$flag=0;
					}
		
	
					if($array[$g]=~/^Sbjct\s+(\d+)\s+([A-Za-z\-]+)\s+(\d+)/)
					{


					push(@subject_seq,$2);
					push(@subject_start,$1);
					push(@subject_end,$3);
					}
		
				}
				@array= ();
				@next=();
				close BLAST1;





		}



	{

		no strict 'refs';
			for (my $i=1;$i<=$#genome_list;$i++)
			{
			my $b=$i-1;
			my $c=$i+1;
			my $var='RESULT$i';
			open($var,"$out_path/process_$species/$name/1_$c")||die("couldnot open the file");
			@{$b}=<$var>;
			close $var;

			}


}


			`mkdir -p $out_path/output_$species/$name`;


						###############writing in common unique############


			open(COMMON,">$out_path/output_$species/$name/total_$name")||die("couldnot open the file");
			open(UNIQUE,">$out_path/output_$species/$name/unique_$name")||die("couldnot open the file");
			open(CONSERVED,">$out_path/output_$species/$name/conserved_with_all_$name")||die("couldnot open the file");


			my $a;


			$a=$genome_list[0];
		
			for(my $r=1;$r<=$#genome_list;$r++)
			{
		

			       $a.=sprintf("\t Hval\t%-15s\t",$genome_list[$r]);
			}

			print CONSERVED "$a\n\n";
			print COMMON "$a\n\n";

			open(DISPENSABLE,">$out_path/output_$species/$name/Dispensable_$name\.txt")||die("couldnot open the file");
			print DISPENSABLE "$a\n\n";
			 
			open(UN,">$out_path/output_$species/$name/unique_genes_$name\.html")||die("could not open file");
			print UN "<HTML><HEAD><h3 align=\"center\">STRAIN-SPECIFIC GENES OF $name</HEAD></h3>\n<h4>LOCATION:PID&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;STRAND&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;GENE_ID&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;SYNONYM&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;COG_CLASS&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ANNOTATION:&nbsp;&nbsp;SEQUENCE</h4><BODY bgcolor= '#E1D5BF' >
<script>
function yuva(a)
{
alert(a);
}


</script>


<pre>\n";
		

			open(AA,"$out_path/process_$species/common/$name.cds")||die("couldnt open file");
			my @genes=<AA>;
			close AA;


	{

		no strict 'refs';

			for(my $i=1;$i<$count_genome;$i++)
			{
				open("VAR$i",">$out_path/output_$species/$name/conserved_in_$i-genomes")||die("couldnot open the file");
				my $v="VAR$i";
				print $v  "$a\n\n";

			}
	}

			open(LIST,"$out_path/process_$species/$name/list")||die("couldnot open the file $!");
			my @list=<LIST>;
			close LIST;	
my ($ref,$ref1);

			for(my $f=0;$f<=$#list;$f++)
			{
				#$hvalue{$h}=();
				@hvalue=();
				#$id{$h}="";
				@geneid=();
				@seq=();
				$flag_val = 0;

				{
					no strict 'refs';
	
					for(my $h=1;$h<=$#genome_list;$h++)
					{
						my $b=$h-1;
					
	
						if(${$b}[$f]=~/^(.+)\s+(\d?\.?\d?\d?)\s+(.+)/)
						{

						$ref=$1;
						$ref1=$ref;
                        $ref=sprintf("%-${max}s",$ref);
						my $h_value=$2;
						my $id_value=$3;
	
						push(@hvalue,$h_value);
						push(@geneid,$id_value);

						}
					}
				}
		      
	
					
					for(my $g=0;$g<=$#hvalue;$g++)
					{
						if($hvalue[$g]>=$hcutoff)
						{
						$flag_val++;
						}
					}
			
					my $hval=join('	',@hvalue);
					my $gid=join('	',@geneid);
	

	
					for(my $s=0;$s<=$#hvalue;$s++)
					{
					my $q=sprintf("% 1.2f",$hvalue[$s]);
					my $w=sprintf("%-15s",$geneid[$s]);

					push(@seq,"\t$q");
					push(@seq,"\t$w");
					}
	

				       if(int($flag_val)==(int($count_genome)-1))
				       {

						print  CONSERVED "$ref @seq\n";
						$count_core++;
						
				       }

			     
			
					if($flag_val>0 && $flag_val<($count_genome-1))
				       {
					print DISPENSABLE " $ref @seq\n";
				       }
	
					{
					no strict 'refs'; no warnings;
					my $filepointer="VAR$flag_val";

		 			print $filepointer "$ref @seq\n";
					}

			    
					print COMMON "$ref\t@seq\n";
	
					
					if($flag_val==0)
					{
						print UNIQUE " $ref\t$hval\n";
						$count_unique++;
						 open(UN,">>$out_path/output_$species/$name/unique_genes_$name\.html")||die("could not open file");
						 
						for(my $x=0;$x<=$#genes;$x++)
						{
						
		
							 if($genes[$x]=~/$ref1/)
							{
						     		chomp ($genes[$x]);
								
								my $ext=$genes[$x+1];
								chomp($ext);
						       		 print UN "$genes[$x]:&nbsp;&nbsp;<a href='' onclick='yuva(\"$ext\")'>sequence</a>\n";

								
							close UN;
		
							  last;
		
							}
					     
						}	
	
				     
					}
					$flag_val="";		
			}


		open(UN,">>$out_path/output_$species/$name/unique_genes_$name\.html")||die("could not open file");
		print UN "</table></BODY>\n</pre>\n</HTML>\n";
		close UN;
		open(OUT,">>$out_path/output_$species/name")||die("couldntopenfile");
		print OUT "$name....conserved in all= $count_core\n";
		print OUT "$name....unique= $count_unique\n";
		print OUT  "...............$count_genome\n";

		$count_unique=0;
		$count_core=0;
		close OUT;
		close UNIQUE;
		close COMMON;
		

		@list=();

}


########################################################################################################33

open(GENOME,"genome_compared")||die("couldnt open file");
@genome_list=<GENOME>;                             #list of the strains in the genome_list array
chomp(@genome_list);
$no=0;

{

no strict 'refs';


		for(my $i=1;$i<=$#genome_list;$i++)
		{
			
			chomp($genome_list[0]);
			chomp($genome_list[$i]);
			my $taken_genome=$genome_list[0];
			my $next_genome=$genome_list[$i];
			my $hash_name="$taken_genome\_$next_genome";
			my $hash_name1="$next_genome\_$taken_genome";

			print  "$hash_name\n";

			
			my $cw=$#genome_list+1+$i;
			open(CORE,">$out_path/process_$species/common/$hash_name\_$hash_name1")||die("couldnt open the file");



			while (my ($key, $value) = each(%$hash_name))
			{
				while (my ($key1, $value1) = each(%$hash_name1))
				{

					if (($key eq $value1) &&($value eq $key1) )
					{
					print CORE "$key....$value\n ";
					push(@{$cw},"$key*$value")
				

					}
				}
			}
		
		
		
		
		}

close CORE;
}

print "core genes identified\n";		
		

$name=$first;	
chomp($name);  
open(TR,"$out_path/process_$species/$name/list")||die ("couldnt open file");  
my @tr=<TR>;
close TR;
open(AA,"$out_path/process_$species/common/$name.cds")||die("couldnt open file");
my @genes=<AA>;
close AA;

			
		
my $cw;	
my $fl;



{
no strict 'refs';


my $a;
		open(FR,">$out_path/process_$species/conseverd")||die ("couldnt open file");  
		for(my $dg=0;$dg<=$#tr;$dg++) 
		{
		$a = "";		
	
		my $fl;
		
		chomp($tr[$dg]);
				for(my $i=1;$i<=$#genome_list;$i++)
				{
					$cw=$#genome_list+1+$i;
					
					for(my $f=1;$f<=$#$cw;$f++)
					{
						if(${$cw}[$f]=~/$tr[$dg](.+)/)
						{
					
						
						my $patt=$1;
						$a=$a.$patt;
						
							
						}
					}
			

				
				}
		print FR "$tr[$dg]*$a\n";
		print "writing conserved gene files\n";
		
		
		}
close FR;	
	
}	

open(FAR,"$out_path/process_$species/conseverd")||die ("couldnt open file");
my @fwr=<FAR>;
close FAR;
open(LIN,">$out_path/output_$species/CORE")||die ("couldnt open file");
my $count_ortholog;

	
	open(CO,">$out_path/output_$species/$name/core_genes_$name\.html")||die("could not open file");
	print CO "<HTML><HEAD><h3 align =\"center\">CORE GENE  OF $name</HEAD></h3>\n<h4>PID&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;PROTIEN_ID&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ANNOTATION:&nbsp;&nbsp;SEQUENCE</h4><BODY bgcolor= '#E1D5BF' >
<script>
function yuva(a)
{
alert(a);
}


</script>
<pre>\n";
			



foreach my $line(@fwr)
{
	my @ae=split('\*',$line);
	my $count=$#ae;
	
	if($line=~/(.+)\*\*/)
	{		
		my $matched_con=$1;	
		if($count==$#genome_list+1)
		{
	
		
				my $ext;
				my $x;
				my @conserved_list;
				for($x=0;$x<=$#genes;$x++)
				{
				
					if($genes[$x]=~/$matched_con/)
					{
						chomp($genes[$x]);						
						$ext=$genes[$x+1];
						chomp($ext);
						last;
					}	
				}
		
				 print CO "$genes[$x]:&nbsp;&nbsp;<a href='' onclick='yuva(\"$ext\")'>sequence</a>\n";	

		





		$line=~s/\*/\t/g;
		
		print LIN "$line";
		push (@conserved_list,$matched_con);
		
		$count_ortholog++;
		}
		
	}
	
	
}






close LIN;
close CO;
my @ae;
 @fwr= ();
 @tr = ();
 @genes =();

print"No of core genes=$count_ortholog\n ";

	 	

			

my (@namevalues,@printcore,@printunique);
@namevalues=@whole;

open(FI,"$out_path/output_$species/name")||die("couldntopenfile");
my @file=<FI>;
	foreach my $e(@whole)

	{
	chomp($e);
	my $matchpat=$e;
		foreach my $l(@file)
		{
			if($l=~/$matchpat\.\.\.\.unique=\s(\d+)/)
		       	{
		 

			my  $uniq=$1;

			   push(@printunique,$uniq);
			}

		}
	}



close FI;
print"process over\n";
@file=();


######################################        LATEX OUTPUT        ############################################



while ($first ne $genome_list[0])
{
$b=shift(@genome_list);                         #shuffling the genome list

push(@genome_list,$b);
}





open (FILEHAN ,">$out_path/output_$species/Result.tex")||die("couldnt open the file");

print FILEHAN "\n\\documentclass{article}
\\usepackage[margin=2.5cm]{geometry}
\\usepackage{tikz,hyperref}
\\title{Core,Dispensable and Strain specific genes plot   }
\\usetikzlibrary{shapes,backgrounds}
\\usepackage{verbatim}
\\usepackage{pgf}
\\hypersetup{
  colorlinks   = true, %Colours links instead of ugly boxes
  urlcolor     = blue, %Colour for external hyperlinks
}
\\date{}
\\begin{document}
\\pagestyle{empty}
\\maketitle
% Now we can draw the sets:
\\begin{tikzpicture}

\\definecolor{orange}{HTML}{FF7F00}\n";
=a

print FILEHAN"\n%\\draw (0,0)--(6,0);
%\\foreach \\x in {0,...,6}
 % \\draw (\\x,0)--(\\x,-.1) node[anchor=east] {\x};

%\\draw (0,0)--(0,6);
%\\foreach \\y in {0,...,6}
% \\draw (0,\\y)--(-.1,\\y) node[anchor=east] {\y};

%\\draw (0,0)--(0,-6);
%\\foreach \\y in {0,...,-6}
 %\\draw (0,\\y)--(-.1,\\y) node[anchor=east] {\y};

%\\draw (0,0)--(-6,0);
%\\foreach \\x in {0,...,-6}
  %\\draw (\\x,0)--(\\x,-.1) node[anchor=east] {\x};\n";
 
=cut
my $xcord=0;
my $ycord=0;
my $angle=360/$number;
$angle=sprintf("%d",$angle);
my $r=3;
my $r1=5;
my $r2=7;
my $r3=8;
my @array=('purple','orange','yellow','red','violet','pink','green','cyan','lime');




for(my $d=0;$d<$number;$d++)
{

	my $color = $array[rand @array];
	my  $theta=($angle*$d*0.01745);
	my  $xneed=$xcord+$r*cos($theta);

	my $yneed=$ycord+$r*sin($theta);
	my $x1=$xcord+$r1*cos($theta);
	my $y1=$ycord+$r1*sin($theta);
	my $x2=$xcord+$r2*cos($theta);
	my $y2=$ycord+$r2*sin($theta);
	

	 $x1=sprintf("%.2f",$x1);
	 $y1=sprintf("%.2f",$y1);
	 $x2=sprintf("%.2f",$x2);
	 $y2=sprintf("%.2f",$y2);
	
	$xneed=sprintf("%.2f",$xneed);
	$yneed=sprintf("%.2f",$yneed);

	my $angl=$angle*$d;
	


	my $printvalue=$printunique[$d];

	my $nameval=$namevalues[$d];

	
	my $namdisplay="../output_$species/$nameval/unique_genes_$nameval\.html";
	chomp($namdisplay);
	my $nameval1=$nameval;
	$nameval1=~s/_/ /;



	print FILEHAN "\n\\draw [fill=$color,rotate around={$angl:($xneed,$yneed)}]($xneed,$yneed) ellipse (3cm and 0.75cm); \n";
	print FILEHAN "\\draw ($x1,$y1)node[]{\\href{$namdisplay}{$printvalue}};\n";
	print FILEHAN "\n\\draw ($x2,$y2)node[]{\\href{http://www.ncbi.nlm.nih.gov/nuccore/$nameval}{$nameval1}};\n";
	
}



my $corevalue=$count_ortholog;
my $nameval=$namevalues[$va];

print FILEHAN "\n\\draw[fill=white] (0,0) circle(1.5);\n";

my $namdis1="../output_$species/$nameval/core_genes_$nameval\.html";
print FILEHAN "\\draw (0,0)node[]{core genes \\href{$namdis1}{$corevalue}};\n";


print FILEHAN "\n\\end{tikzpicture}

\n\n\n\n\n";


$nameval=$namevalues[$va];
print FILEHAN"\n\\begin{flushright}{\\href{../output_$species/$nameval/Dispensable_$nameval\.txt}{Click here for Dispensable gene}}\\end{flushright}\n";
print FILEHAN"\n\\begin{flushright}{\\href{../output_$species/no_of_genes.txt}{Number of genes in each genome}}\\end{flushright}\n";




print FILEHAN"\\null
\\vfill


*  click on the numbers to see the unique gene sequences\n
*  click on the NCBI-IDs to go to NCBI genome page\n


\n\\end{document}";
close FILEHAN;

`pdflatex --output-directory=$out_path/output_$species/ $out_path/output_$species/Result.tex`;

my $error = $?;

if($error == 0)
{
	print "LateX pdf succesfuly created\n";
}
else
{
	print "Latex file not created";
}





			
}
