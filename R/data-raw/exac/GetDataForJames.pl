use strict;
use warnings;

use POSIX;
use DBI;

#--> Connect to database
my $dbm = DBI->connect('DBI:mysql:ICC_MUTATIONS:cmr000236', 'icc', 'A024N') or die ("cant connect $DBI::errstr\n");

# Get mutation table where ExAC count >0 with ClinVar mset id if one exists for the 8 sarcomeric genes
print "Looking at mutation table\n";

$"="\t";
open(MUT, ">Mutation_inExAC_sarcomericGenesOnly.txt");

my $col_count=0;
my $sql0=$dbm->prepare("SHOW COLUMNS FROM mutation"); 
$sql0->execute  or die ("cant execute statement $DBI::errstr\n");
while (my $arrayref0 = $sql0->fetchrow_arrayref)
{
	my $colname=$arrayref0->[0];
	if ($col_count==0)
	{
		print MUT "$colname";
	}
	else
	{
		print MUT "\t$colname";
	}
	$col_count++;
}
print MUT "\tClinVar_mset\tClinVar_sig\n";

my $sql1=$dbm->prepare("SELECT * FROM mutation WHERE mut_exac_count>0 AND mut_gene_id IN ('25','30','38','41','26','29','31','32')"); 
$sql1->execute  or die ("cant execute statement $DBI::errstr\n");
while (my $arrayref = $sql1->fetchrow_arrayref)
{
	my $length = scalar @{ $arrayref };
	my $mset="";
	my $class="";
	my $mut_id="";
	for (my $n=0; $n<$length; $n++)
	{
		my $data;
		$data=$arrayref->[$n];
		if ($n==0)
		{
			$mut_id=$data;
		}
		else
		{
			print MUT "\t";
		}
		if ($data){}
		else
		{
			$data="";
		}
		print MUT "$data";
	}
	my $count=0;
	my $sql2=$dbm->prepare("SELECT mc_mset_id, cr_clin_sig FROM mut_clinvar, clinvar_reports, disease_class WHERE mc_mut_id='$mut_id' AND mc_mset_id=cr_mset_id AND cr_trait=dc_description AND dc_subclass='HCM'"); 
	$sql2->execute  or die ("cant execute statement $DBI::errstr\n");
	while (my $arrayref2 = $sql2->fetchrow_arrayref)
	{
		$mset=$arrayref2->[0];
		my $clin_sig=$arrayref2->[1];
		if ($count==0)
		{
			$class=$clin_sig;
		}
		else
		{
			$class=$class.";".$clin_sig;
		}
		$count++;
	}
	print MUT "\t$mset\t$class\n";
}

# Get the gene table for the 8 sarcomeric genes

print "Looking at gene table\n";

open(GEN, ">GeneInfo.txt");

my $col_count2=0;
my $sql3=$dbm->prepare("SHOW COLUMNS FROM gene"); 
$sql3->execute  or die ("cant execute statement $DBI::errstr\n");
while (my $arrayref3 = $sql3->fetchrow_arrayref)
{
	my $colname=$arrayref3->[0];
	if ($col_count2==0)
	{
		print GEN "$colname";
	}
	else
	{
		print GEN "\t$colname";
	}
	$col_count2++;
}
print GEN "\n";

my $sql2=$dbm->prepare("SELECT * FROM gene WHERE gene_id IN ('25','30','38','41','26','29','31','32')"); 
$sql2->execute  or die ("cant execute statement $DBI::errstr\n");
while (my $arrayref2 = $sql2->fetchrow_arrayref)
{
	my $length2 = scalar @{ $arrayref2 };
	for (my $n=0; $n<$length2; $n++)
	{
		my $data;
		$data=$arrayref2->[$n];
		if ($n==0){}
		else
		{
			print GEN "\t";
		}
		if ($data){}
		else
		{
			$data="";
		}
		print GEN "$data";
	}
	print GEN "\n";
}

# Get the ExAC data and population tables

print "Looking at ExAC data\n";

open(EXP, ">ExACPopulations.txt");
print EXP "ep_id\tep_code\tep_name\tep_max_alleles\n";

my $sql4=$dbm->prepare("SELECT ep_id, ep_code, ep_name, ep_max_alleles FROM exac_population"); 
$sql4->execute  or die ("cant execute statement $DBI::errstr\n");
while (my $arrayref4 = $sql4->fetchrow_arrayref)
{
	my $ep_id=$arrayref4->[0];
	my $ep_code=$arrayref4->[1];
	my $ep_name=$arrayref4->[2];
	my $ep_max_alleles=$arrayref4->[3];

	print EXP "$ep_id\t$ep_code\t$ep_name\t$ep_max_alleles\n";
}

open (EXD, ">ExACPASSData.txt");

my $col_count3=0;
my $sql5=$dbm->prepare("SHOW COLUMNS FROM exac_data"); 
$sql5->execute  or die ("cant execute statement $DBI::errstr\n");
while (my $arrayref5 = $sql5->fetchrow_arrayref)
{
	my $colname=$arrayref5->[0];
	if ($col_count3==0)
	{
		print EXD "$colname";
	}
	else
	{
		print EXD "\t$colname";
	}
	$col_count3++;
}
print EXD "\n";

my $sql6=$dbm->prepare("SELECT * FROM exac_data"); 
$sql6->execute  or die ("cant execute statement $DBI::errstr\n");
while (my $arrayref6 = $sql6->fetchrow_arrayref)
{
	my $length3 = scalar @{ $arrayref6 };
	for (my $n=0; $n<$length3; $n++)
	{
		my $data;
		$data=$arrayref6->[$n];
		if ($n==0){}
		else
		{
			print EXD "\t";
		}
		if ($data){}
		else
		{
			$data="0";
		}
		print EXD "$data";
	}
	print EXD "\n";
}


