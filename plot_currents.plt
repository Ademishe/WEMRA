#!/usr/local/bin/gnuplot
set title "8 периодов РДГ"
set terminal pdfcairo enhanced
set output "RDG1.pdf"
set encoding utf8
set xlabel "Время" font "Helvetica,15"
set ylabel "Выходная мощность" font "Helvetica,15"
set grid xtics ytics

plot "power.dat" u 8:11 w linespoints ps 0.2
