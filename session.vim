let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/uni/fys4411/project1/src
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +1 Hamiltonians/hamiltonian.h
badd +0 Hamiltonians/harmonicoscillator.h
badd +0 InitialStates/initialstate.h
badd +0 InitialStates/randomuniform.h
badd +0 Math/random.h
badd +0 particle.h
badd +0 sampler.h
badd +0 system.h
badd +0 WaveFunctions/simplegaussian.h
badd +0 WaveFunctions/wavefunction.h
badd +48 Hamiltonians/hamiltonian.cpp
badd +46 Hamiltonians/harmonicoscillator.cpp
badd +0 InitialStates/initialstate.cpp
badd +0 InitialStates/randomuniform.cpp
badd +1 main.cpp
badd +0 particle.cpp
badd +88 sampler.cpp
badd +0 system.cpp
badd +0 WaveFunctions/simplegaussian.cpp
badd +0 WaveFunctions/wavefunction.cpp
argglobal
%argdel
$argadd Hamiltonians/hamiltonian.h
$argadd Hamiltonians/harmonicoscillator.h
$argadd InitialStates/initialstate.h
$argadd InitialStates/randomuniform.h
$argadd Math/random.h
$argadd particle.h
$argadd sampler.h
$argadd system.h
$argadd WaveFunctions/simplegaussian.h
$argadd WaveFunctions/wavefunction.h
$argadd build/CMakeFiles/3.19.4/CompilerIdCXX/CMakeCXXCompilerId.cpp
$argadd Hamiltonians/hamiltonian.cpp
$argadd Hamiltonians/harmonicoscillator.cpp
$argadd InitialStates/initialstate.cpp
$argadd InitialStates/randomuniform.cpp
$argadd main.cpp
$argadd particle.cpp
$argadd sampler.cpp
$argadd system.cpp
$argadd WaveFunctions/simplegaussian.cpp
$argadd WaveFunctions/wavefunction.cpp
edit main.cpp
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
exe 'vert 1resize ' . ((&columns * 74 + 85) / 170)
exe 'vert 2resize ' . ((&columns * 95 + 85) / 170)
argglobal
if bufexists("main.cpp") | buffer main.cpp | else | edit main.cpp | endif
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let s:l = 52 - ((26 * winheight(0) + 22) / 45)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
52
normal! 096|
wincmd w
argglobal
if bufexists("Hamiltonians/hamiltonian.cpp") | buffer Hamiltonians/hamiltonian.cpp | else | edit Hamiltonians/hamiltonian.cpp | endif
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let s:l = 48 - ((40 * winheight(0) + 22) / 45)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
48
normal! 019|
wincmd w
2wincmd w
exe 'vert 1resize ' . ((&columns * 74 + 85) / 170)
exe 'vert 2resize ' . ((&columns * 95 + 85) / 170)
tabnext 1
if exists('s:wipebuf') && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 winminheight=1 winminwidth=1 shortmess=filnxtToOF
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
