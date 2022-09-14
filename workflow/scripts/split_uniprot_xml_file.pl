#!/usr/bin/env perl

use strict;
use warnings;
use Cwd qw(cwd);
use autodie qw(:all);
use File::Basename qw(fileparse);
use File::Path qw(make_path);
use File::Spec;
use Getopt::Long qw(:config auto_help auto_version);
use PerlIO::gzip;
use Pod::Usage qw(pod2usage);
use XML::LibXML::Reader;

# unbuffer error and output streams (make sure STDOUT
# is last so that it remains the default filehandle)
select(STDERR); $| = 1;
select(STDOUT); $| = 1;

my %xml_type_name = (
    &XML_READER_TYPE_ELEMENT     => 'ELEMENT',
    &XML_READER_TYPE_END_ELEMENT => 'END_ELEMENT',
);

my $xmlns = "http://uniprot.org/uniprot";

my $xml_header = <<"XML_HEADER";
<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="$xmlns"
 xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
 xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/docs/uniprot.xsd">
XML_HEADER

my $xml_footer = <<"XML_FOOTER";
<copyright>
Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms Distributed under the Creative Commons Attribution (CC BY 4.0) License
</copyright>
</uniprot>
XML_FOOTER

sub write_split_file {
    my ($elems_ref, $out_dir, $file_basename, $num, $sep) = @_;
    $num = sprintf("%04d", $num);
    make_path($out_dir, { chmod => 0755 }) if !-e $out_dir;
    my $split_filename = "${file_basename}_${num}.xml.gz";
    print "Writing $split_filename\n";
    open(my $fh, '>:gzip', File::Spec->catfile($out_dir, $split_filename));
    print $fh $xml_header;
    print $fh join($sep, @{$elems_ref});
    print $fh $sep . $xml_footer;
    close($fh);
}

my $in_file = '';
my $out_dir = cwd();
my $file_basename = 'uniprot';
my $split_size = 1000000;
my $parser = 'perl_regex';
my $verbose = 0;
GetOptions(
    'in-file|i:s' => \$in_file,
    'out-dir|0:s' => \$out_dir,
    'basename|b:s' => \$file_basename,
    'split-size|s:i' => \$split_size,
    'parser|p:s' => \$parser,
) || pod2usage(-verbose => 0);
if ($in_file eq '') {
    pod2usage(-message => "Required --in-file");
}
elsif (!-f $in_file) {
    pod2usage(-message => "Invalid --in-file: $in_file");
}

my @split_elems;
my $split_num = 1;
my $num_split_elems = 0;
my $num_total_elems = 0;
my $elem_sep = $parser eq "perl_lxml" ? "" : "\n";
open(my $xml_fh, '<:gzip', $in_file);
if ($parser eq "perl_lxml") {
    my $reader = XML::LibXML::Reader->new(IO => $xml_fh);
    while($reader->read) {
        next unless $reader->nodeType == XML_READER_TYPE_ELEMENT;
        next unless $reader->name eq 'entry';
        push(@split_elems, $reader->readOuterXml);
        $num_split_elems++;
        $num_total_elems++;
        if ($num_split_elems == $split_size) {
            write_split_file(\@split_elems, $out_dir, $file_basename, $split_num, $elem_sep);
            $split_num++;
            @split_elems = ();
            $num_split_elems = 0;
        }
        $reader->next;
    }
    $reader->close;
}
else {
    my $in_elem = 0;
    my @elem_lines;
    my $elem_s_regex = qr/\s*<\s*entry\s+.*?>/;
    my $elem_e_regex = qr/\s*<\/\s*entry\s*>/;
    while (<$xml_fh>) {
        s/\s+$//o;
        if (/$elem_s_regex/) {
            $in_elem = 1;
            push(@elem_lines, $_);
        }
        elsif (/$elem_e_regex/) {
            die "Closing entry tag without previously reading opening tag" if !$in_elem;
            push(@elem_lines, $_);
            $in_elem = 0;
            push(@split_elems, join("\n", @elem_lines));
            @elem_lines = ();
            $num_split_elems++;
            $num_total_elems++;
            if ($num_split_elems == $split_size) {
                write_split_file(\@split_elems, $out_dir, $file_basename, $split_num, $elem_sep);
                $split_num++;
                @split_elems = ();
                $num_split_elems = 0;
            }
        }
        elsif ($in_elem) {
            push(@elem_lines, $_);
        }
    }
}
close($xml_fh);
write_split_file(\@split_elems, $out_dir, $file_basename, $split_num, $elem_sep);
print "Parsed $num_total_elems $file_basename records\n";
exit(0);

__END__

=head1 NAME

split_uniprot_xml_file.pl - Split UniProt XML file into parts

=head1 SYNOPSIS

 split_uniprot_xml_file.pl [options]

 Options:
    --in-file|-i <file>     Input UniProt gzipped XML file
                            (required)
    --out-dir|-o <dir>      Output directory
                            (default = current working directory)
    --basename|-b <str>     Output file basename
                            (required)
    --split-size|-s <int>   Number of elements in each split file
                            (default = 1000000)
    --parser|-p <str>       Parser type
                            (default = perl_regex)
    --help                  Display usage and exit
    --version               Display program version and exit

=cut