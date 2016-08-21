#!/usr/bin/perl
# 
# Take text from argv[1] and convert to Ｆｕｌｌｗｉｄｔｈ.
# It's not all that hard to do, it's just remapping ASCII to unicode points.
# Using a website for something this simple is for n00bs.
# This text attempts to print UTF-8 to your console, or to a file if you ask
# your shell to redirect to a file.

use utf8;
use encoding 'UTF-8';

my $ascii = $ARGV[0];
die unless defined($ascii);
die if $ascii eq '';

binmode(STDOUT,":utf8");

sub makewide($) {
    my $a = $_[0];

       if ($a eq '!') { return '！'; }
    elsif ($a eq '\"') { return '＂'; }
    elsif ($a eq '#') { return '＃'; }
    elsif ($a eq '$') { return '＄'; }
    elsif ($a eq '%') { return '％'; }
    elsif ($a eq '&') { return '＆'; }
    elsif ($a eq '\'') { return '＇'; }
    elsif ($a eq '(') { return '（'; }
    elsif ($a eq ')') { return '）'; }
    elsif ($a eq '*') { return '＊'; }
    elsif ($a eq '+') { return '＋'; }
    elsif ($a eq ',') { return '，'; }
    elsif ($a eq '-') { return '－'; }
    elsif ($a eq '.') { return '．'; }
    elsif ($a eq '/') { return '／'; }

    elsif ($a eq ':') { return '：'; }
    elsif ($a eq ';') { return '；'; }
    elsif ($a eq '<') { return '＜'; }
    elsif ($a eq '=') { return '＝'; }
    elsif ($a eq '>') { return '＞'; }
    elsif ($a eq '?') { return '？'; }

    elsif ($a eq '@') { return '＠'; }

    elsif ($a eq '[') { return '［'; }
    elsif ($a eq '\\') { return '＼'; }
    elsif ($a eq ']') { return '］'; }
    elsif ($a eq '^') { return '＾'; }
    elsif ($a eq '_') { return '＿'; }

    elsif ($a eq '`') { return '｀'; }
    elsif ($a eq '{') { return '｛'; }
    elsif ($a eq '|') { return '｜'; }
    elsif ($a eq '}') { return '｝'; }
    elsif ($a eq '~') { return '～'; }

    elsif ($a eq '.') { return '｡'; }

    elsif ($a eq ' ') { return '　'; }

    $ac = ord($a);
    if ($ac >= ord('A') && $ac <= ord('Z')) { return chr($ac + ord('Ａ') - ord('A')); }
    if ($ac >= ord('a') && $ac <= ord('z')) { return chr($ac + ord('ａ') - ord('a')); }
    if ($ac >= ord('0') && $ac <= ord('9')) { return chr($ac + ord('０') - ord('0')); }

    return $a;
}

$final = '';
for ($i=0;$i < length($ascii);$i++) {
    $ci = substr($ascii,$i,1);
    $co = makewide($ci);
    $final .= $co;
}

print "$final\n";

