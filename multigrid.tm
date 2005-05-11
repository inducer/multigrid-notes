<TeXmacs|1.0.5>

<style|<tuple|book|number-long-article>>

<\body>
  <with|mode|math|<assign|ip|<macro|left|right|(<arg|left>,<arg|right>)>><assign|l2norm|<macro|arg|<left|\|\|><arg|arg><right|\|\|>>><assign|enorm|<macro|arg|<left|\|><left|\|\|><arg|arg><right|\|\|><right|\|>>>>Multigrid
  Methods

  Prof. Christian Wieners, SS 2005

  Korrekturen und Vorschläge an <with|font-family|tt|kloeckner@math.uni-karlsruhe.de>.

  <\table-of-contents|toc>
    <vspace*|1fn><with|font-series|bold|math-font-series|bold|Table of
    contents> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-1><vspace|0.5fn>

    1<space|2spc>Introduction <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-2>

    <with|par-left|1.5fn|1.1<space|2spc>The model problem
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-3>>

    <with|par-left|1.5fn|1.2<space|2spc>Iterative Methods for Linear Systems
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-4>>

    <with|par-left|1.5fn|1.3<space|2spc>The Multigrid Idea
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-5>>

    2<space|2spc>A two-level method <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-6>

    <with|par-left|1.5fn|2.1<space|2spc>The model problem
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-7>>

    <with|par-left|1.5fn|2.2<space|2spc>Implementation of the two-level
    method <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-11>>

    3<space|2spc>Classical two-level analysis
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-12>
  </table-of-contents>

  <section|Introduction>

  <subsection|The model problem>

  <label|subsec:fd-model-problem>On <with|mode|math|\<Omega\>=(0,1)<rsup|2>>,
  consider the Laplace problem.

  <\quotation>
    Find <with|mode|math|u:<wide|\<Omega\>|\<bar\>>\<rightarrow\>\<bbb-R\>>
    such that

    <\eqnarray*>
      <tformat|<table|<row|<cell|-\<Delta\>u>|<cell|=>|<cell|f<space|1fn>(x\<in\>\<Omega\>)>>|<row|<cell|u>|<cell|=>|<cell|0<space|1fn>(x\<in\>\<partial\>\<Omega\>)>>>>
    </eqnarray*>
  </quotation>

  For <with|mode|math|n\<in\>\<bbb-N\>>, <with|mode|math|h=1/n> and
  <with|mode|math|\<cal-N\><rsub|h>=h\<bbb-Z\><rsup|2>\<cap\>\<Omega\>>. We
  set

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|1,h><rsup|+>u(x<rsub|1>,x<rsub|2>)>|<cell|\<assign\>>|<cell|<frac|1|h>(u(x<rsub|1>+h,x<rsub|2>)-u(x<rsub|1>,x<rsub|2>))>>|<row|<cell|\<partial\><rsub|1,h><rsup|->u(x<rsub|1>,x<rsub|2>)>|<cell|\<assign\>>|<cell|<frac|1|h>(u(x<rsub|1>-h,x<rsub|2>)-u(x<rsub|1>,x<rsub|2>))>>|<row|<cell|\<partial\><rsub|1,h><rsup|2>u(x<rsub|1>,x<rsub|2>)>|<cell|\<assign\>>|<cell|<frac|1|h<rsup|2>>(u(x<rsub|1>-h,x<rsub|2>)-2u(x<rsub|1>,x<rsub|2>)+u(x<rsub|1>+h,x<rsub|2>))>>|<row|<cell|-\<Delta\><rsub|h>u(x<rsub|1>,x<rsub|2>)>|<cell|\<assign\>>|<cell|<frac|1|h<rsup|2>>(4u(x<rsub|1>,x<rsub|2>)-u(x<rsub|1>-h,x<rsub|2>)-u(x<rsub|1>+h,x<rsub|2>)-u(x<rsub|1>,x<rsub|2>+h)-u(x<rsub|1>,x<rsub|2>-h))>>|<row|<cell|-\<Delta\><rsub|h>>|<cell|=>|<cell|<left|[><tabular*|<tformat|<table|<row|<cell|>|<cell|-1>|<cell|>>|<row|<cell|-1>|<cell|4>|<cell|-1>>|<row|<cell|>|<cell|-1>|<cell|>>>>><right|]><space|1fn><with|mode|text|(Finite
    difference star)>>>>>
  </eqnarray*>

  Find a grid function <with|mode|math|u:<wide|\<cal-N\><rsub|h>|\<bar\>>\<assign\>\<bbb-Z\><rsup|2>\<cap\><wide|\<Omega\>|\<bar\>>\<rightarrow\>\<bbb-R\>>
  s.t. <with|mode|math|-\<Delta\><rsub|h>u<rsub|h>=f<rsub|h>> in
  <with|mode|math|\<cal-N\><rsub|h>>, <with|mode|math|u<rsub|h>=0> on
  <with|mode|math|\<partial\>\<cal-N\><rsub|h>=<wide|\<cal-N\><rsub|h>|\<bar\>>\<setminus\>\<cal-N\><rsub|h>>.
  By lexicographic ordering on <with|mode|math|\<cal-N\><rsub|h>={(i*h,j*h):i,j=1,\<ldots\>,n-1}>,
  <with|mode|math|u<rsub|h>\<in\>\<cal-N\><rsub|h>> is represented as
  <with|mode|math|\<b-u\>\<in\>\<bbb-R\><rsup|n>>
  (<with|mode|math|N=(n-1)<rsup|2>)>: <with|mode|math|\<b-u\>[i(n-1)+j]=u<rsub|h>(i*h,j*h)>.

  <\equation*>
    \<b-A\>\<b-u\>=\<b-f\>, <with|mode|text|with>
    \<b-A\>=<matrix|<tformat|<table|<row|<cell|T>|<cell|-I>|<cell|>|<cell|>|<cell|>>|<row|<cell|-I>|<cell|T>|<cell|-I>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|\<ddots\>>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|-I>|<cell|T>|<cell|-I>>|<row|<cell|>|<cell|>|<cell|>|<cell|-I>|<cell|T>>>>>\<in\>\<bbb-R\><rsup|N,N>\<backsimeq\>(\<bbb-R\><rsup|n-1,n-1>)<rsup|n-1,n-1>,
  </equation*>

  <\equation*>
    T=<matrix|<tformat|<table|<row|<cell|4>|<cell|-1>|<cell|>|<cell|>|<cell|>>|<row|<cell|-1>|<cell|4>|<cell|-1>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|\<ddots\>>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|-1>|<cell|4>|<cell|-1>>|<row|<cell|>|<cell|>|<cell|>|<cell|-1>|<cell|4>>>>>\<in\>\<bbb-R\><rsup|n-1,n-1>,
  </equation*>

  <with|mode|math|I=diag(1,\<ldots\>,1)\<in\>\<bbb-R\><rsup|n-1,n-1>>.

  <\definition>
    <with|mode|math|\<b-A\>> sparse<with|mode|math|:\<Leftrightarrow\>><with|mode|math|\<b-A\>\<b-u\>>
    requires <with|mode|math|O(n)> operations.
  </definition>

  Here: max. 5 entries per row, bandwidth
  <with|mode|math|n><with|mode|math|\<Rightarrow\>><with|mode|math|\<b-A\>\<b-u\>>
  requires <with|mode|math|5n> operations. <with|mode|math|\<b-A\><rsup|-1>>
  requires <with|mode|math|O(N*n<rsup|2>)=O(N<rsup|2>)=O(n<rsup|4>)>
  operations and <with|mode|math|O(N*n)> storage (worse in 3D!).

  <subsection|Iterative Methods for Linear Systems>

  Let <with|mode|math|\<b-A\>\<in\>\<bbb-R\><rsup|N,N>> be regular,
  <with|mode|math|\<b-f\>\<in\>\<bbb-R\><rsup|N>> the right-hand side, and
  <with|mode|math|\<b-B\>\<in\>\<bbb-R\><rsup|N,N>> regular (the
  ``preconditioner''). For a given <with|mode|math|\<b-u\><rsup|0>\<in\>\<bbb-R\><rsup|N>>
  consider for <with|mode|math|k=1,2,\<ldots\>>
  <with|mode|math|\<b-u\><rsup|k>=\<b-u\><rsup|k-1>+\<b-B\>(\<b-f\>-\<b-A\>\<b-u\><rsup|k-1>)>.
  The error is, using <with|mode|math|\<b-u\>=\<b-A\><rsup|-1>\<b-f\>>,

  <\equation*>
    \<b-u\><rsup|k>-\<b-u\>=\<b-u\><rsup|k-1>-\<b-u\>+\<b-B\>(\<b-A\>\<b-u\>-\<b-A\>\<b-u\><rsup|k-1>)=(\<b-I\>-\<b-B\>\<b-A\>)(\<b-u\><rsup|k-1>-\<b-u\>)=(\<b-I\>-\<b-B\>\<b-A\>)<rsup|k>(\<b-u\><rsup|0>-\<b-u\>)
  </equation*>

  <\definition>
    <with|mode|math|\<sigma\>(\<b-A\>)\<assign\>{\<lambda\>\<in\>\<bbb-C\>:det(\<b-A\>-\<lambda\>\<b-I\>)=0}>,
    called the spectrum. <with|mode|math|\<rho\>(\<b-A\>)\<assign\>max{\|\<lambda\>\|:\<lambda\>\<in\>\<sigma\>(\<b-A\>)}>,
    called the spectral radius. <with|mode|math|\<kappa\>(\<b-A\>)\<assign\>\<rho\>(\<b-A\>)\<rho\>(\<b-A\><rsup|-1>)>,
    called the condition number. <with|mode|math|\<\|\|\>\<b-A\>\<\|\|\>\<assign\>sup<rsub|\<\|\|\>\<b-u\>\<\|\|\>=1>\<\|\|\>\<b-A\>\<b-u\>\<\|\|\>>,
    where <with|mode|math|\<\|\|\>\<b-u\>\<\|\|\>\<assign\><sqrt|u<rsup|T>u>>,
    called the spectral norm.
  </definition>

  <\theorem>
    <\equation*>
      lim<rsub|k\<rightarrow\>\<infty\>>\<\|\|\>(\<b-I\>-\<b-B\>\<b-A\>)<rsup|k>\<b-e\>\<\|\|\>=0
    </equation*>

    for any <with|mode|math|\<b-e\>\<in\>\<bbb-R\><rsup|N>><with|mode|math|\<Leftrightarrow\>><with|mode|math|\<rho\>(\<b-I\>-\<b-B\>\<b-A\>)\<less\>1>.
  </theorem>

  Application to <with|mode|math|\<b-B\>=\<theta\>\<b-I\>> with
  <with|mode|math|0\<less\>\<theta\>\<less\>1/\<rho\>(\<b-A\>)>:
  <with|mode|math|\<sigma\>(\<b-I\>-\<theta\>\<b-A\>)\<subset\>[-\<theta\>\<sigma\>(\<b-A\>),\<theta\>\<sigma\>(\<b-A\>)]><with|mode|math|\<Rightarrow\>>Iteration
  converges.

  <subsection|The Multigrid Idea>

  The basic assumption is that the fine-grid iteration has a hard time
  eliminating low-frequency errors. We consider the model problem. Define
  <with|mode|math|w<rsub|k><rsup|i>=<sqrt|h/2>*sin(h*k*i\<pi\>)>
  (<with|mode|math|i,k=1,\<ldots\>,n-1>, <with|mode|math|h=1/n>). Next,
  define <with|mode|math|\<b-w\><rsup|i,j>\<in\>\<bbb-R\><rsup|N>> by
  <with|mode|math|\<b-w\><rsup|i,j>[k(n-1)+l]=w<rsup|i><rsub|k>w<rsup|j><rsub|l>>.

  <\lemma>
    <\itemize>
      <item><with|mode|math|\<b-A\>\<b-w\><rsup|i,j>=\<lambda\><rsup|i,j>\<b-w\><rsup|i,j>>
      with <with|mode|math|\<lambda\><rsup|i,j>=4h<rsup|-2>(sin(i\<pi\>h/2)+sin<rsup|2>(j\<pi\>h/2))>,

      <item><with|mode|math|(\<b-w\><rsup|i,j>)> form a complete ONB,

      <item><with|mode|math|\<rho\>(\<b-A\>)=\<lambda\><rsub|max>=\<lambda\><rsup|n-1,n-1>=8h<rsup|-2>cos(\<pi\>h/2)=O(h<rsup|-2>)>,

      <with|mode|math|\<lambda\><rsub|min>=\<lambda\><rsup|1,1>=8h<rsup|-2>sin<rsup|2>(\<pi\>h/2)=O(1)>.
    </itemize>
  </lemma>

  Application to <with|mode|math|\<b-I\>-\<theta\>\<b-A\>>:
  <with|mode|math|\<sigma\>(\<b-I\>-\<theta\>\<b-A\>)\<subset\>[1-\<theta\>\<lambda\><rsub|max>,1-\<theta\>\<lambda\><rsub|min>]>,
  \ <with|mode|math|\<rho\>(\<b-I\>-\<theta\>\<b-A\>)=max{\|1-\<theta\>\<lambda\><rsub|max>\|,\|1-\<theta\>\<lambda\><rsub|min>\|}>.
  The optimal choice is

  <\equation*>
    \<theta\><rsup|\<ast\>>=<frac|2|\<lambda\><rsub|max>+\<lambda\><rsub|min>>=<frac|2|8h<rsup|-2>>=<frac|1|4>h<rsup|2>.
  </equation*>

  <with|mode|math|\<Rightarrow\>\<rho\>(\<b-I\>-\<theta\>\<b-A\>)=1-O(h<rsup|-2>)><with|mode|math|\<rightarrow\>>slow!

  <\lemma>
    Let <with|mode|math|\<cal-F\>={(i,j):max{i,j}\<geqslant\>n/2}>,
    <with|mode|math|X=span{\<b-w\><rsup|i,j>:(i,j)\<in\>\<cal-F\>}>, called
    the subspace of strongly oscillating ``functions'', <with|mode|math|dim
    X=3N/4>. Set <with|mode|math|\<theta\>=h<rsup|2>/8>, and assume for the
    initial error <with|mode|math|\<b-e\><rsup|0>=\<b-u\><rsup|0>-\<b-u\>\<in\>X>.
    Then, we have <with|mode|math|\<b-e\><rsup|k>=(\<b-I\>-\<theta\>\<b-A\>)\<b-e\><rsup|k-1>\<in\>X>
    and <with|mode|math|\<\|\|\>\<b-e\><rsup|k>\<\|\|\>\<leqslant\>3/4\<\|\|\>\<b-e\><rsup|k-1>\<\|\|\>\<leqslant\>(3/4)<rsup|k>\<\|\|\>\<b-e\><rsup|0>\<\|\|\>>.
  </lemma>

  <\proof>
    Consider

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<\|\|\>\<b-e\><rsup|k>\<\|\|\><rsup|2>>|<cell|=>|<cell|<big|sum><rsub|i,j>(\<b-e\><rsup|k>\<cdot\>\<b-w\><rsup|i,j>)<rsup|2>=<big|sum><rsub|i,j><left|[>(\<b-I\>-\<theta\>\<b-A\>)<rsup|k>\<b-e\><rsup|0>\<cdot\>\<b-w\><rsup|i,j><right|]><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|i,j><left|[>\<b-e\><rsup|0>\<cdot\>(\<b-I\>-\<theta\>\<b-A\>)<rsup|k>\<b-w\><rsup|i,j><right|]><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|(i,j)\<nin\>\<cal-F\>>(1-\<theta\>\<lambda\><rsup|i,j>)<rsup|2k><wide*|<left|[><wide*|\<b-e\><rsup|0>\<cdot\>\<b-w\><rsup|i,j>|\<wide-underbrace\>><rsub|\<b-e\><rsup|0>\<in\>X\<perp\>{\<b-w\><rsup|i,j>:i,j\<nin\>\<cal-F\>}
      ><right|]><rsup|2>|\<wide-underbrace\>><rsub|=0>+<big|sum><rsub|(i,j)\<in\>\<cal-F\>>(1-\<theta\>\<lambda\><rsup|i,j>)<rsup|2k><left|[>\<b-e\><rsup|0>\<cdot\>\<b-w\><rsup|i,j><right|]><rsup|2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<left|(>max<rsub|(i,j)\<in\>\<cal-F\>>(1-\<theta\>\<lambda\><rsup|i,j>)<rsup|2><right|)><rsup|k>\<\|\|\>\<b-e\><rsup|0>\<\|\|\>.>>>>
    </eqnarray*>

    <\eqnarray*>
      <tformat|<table|<row|<cell|{\<lambda\><rsup|i,j>:(i,j)\<in\>\<cal-F\>}>|<cell|=>|<cell|4h<rsup|-2>{sin<rsup|2>(\<pi\>i*h/2)+sin<rsup|2>(\<pi\>j*h/2):i=1,\<ldots\>,n-1,j=n/2,\<ldots\>,n-1}>>|<row|<cell|>|<cell|\<subset\>>|<cell|4h<rsup|-2><left|[>sin<rsup|2>(\<pi\>/4),2]=h<rsup|-2>[2,8]>>>>
    </eqnarray*>

    <\eqnarray*>
      <tformat|<table|<row|<cell|max<rsub|(i,j)\<in\>\<cal-F\>>\|1-\<theta\>\<lambda\><rsup|i,j>\|>|<cell|=>|<cell|\|1-<frac|1|8>h<rsup|2>\<lambda\><rsup|1,n/2>\|>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|max<rsub|t\<in\>h<rsup|-2>[2,8]>\|1-<frac|1|8>h<rsup|2>t\|=<frac|3|4>>>>>
    </eqnarray*>
  </proof>

  Consequence: Relaxation <with|mode|math|\<b-I\>-\<theta\>\<b-A\>>
  ``smoothing'' the error, slow for the smooth error components. Problem:
  <with|mode|math|X\<subset\>\<bbb-R\><rsup|N>> is not known in general.
  Idea: Continue relaxation with <em|coarse grid correction>.

  <section|A two-level method>

  <subsection|The model problem>

  We consider a weak formulation of the Laplace problem

  <\eqnarray*>
    <tformat|<table|<row|<cell|-\<Delta\>u>|<cell|=>|<cell|f<space|1fn><with|mode|text|on
    <with|mode|math|\<Omega\>\<subset\>\<bbb-R\><rsup|2>>>>>|<row|<cell|u>|<cell|=>|<cell|0<space|1fn><with|mode|text|on
    <with|mode|math|\<partial\>\<Omega\>>>>>>>
  </eqnarray*>

  Let <with|mode|math|v:\<Omega\>\<rightarrow\>\<bbb-R\>> be a test function
  with <with|mode|math|v\|<rsub|\<partial\>\<Omega\>>=0>. Gauÿ:
  <with|mode|math|\<sigma\>\<assign\>v\<nabla\>u>

  <\equation*>
    <big|int><rsub|\<Omega\>>div \<sigma\>*d*v=<big|int><rsub|\<partial\>\<Omega\>>\<sigma\>\<cdot\>n*d*a=0
  </equation*>

  <\equation*>
    div \<sigma\>=\<nabla\>v\<cdot\>\<nabla\>u+v\<Delta\>u=\<nabla\>u\<cdot\>\<nabla\>v-f*v
  </equation*>

  <\equation*>
    \<Rightarrow\><big|int><rsub|\<Omega\>>\<nabla\>u\<cdot\>\<nabla\>v=<big|int>f*v*d*v
  </equation*>

  <\definition>
    Let

    <\eqnarray*>
      <tformat|<table|<row|<cell|<ip|v|w>>|<cell|\<assign\>>|<cell|<big|int><rsub|\<Omega\>>v*w*d*v<space|1fn><with|mode|text|(inner
      product in <with|mode|math|L<rsup|2>(\<Omega\>)>)>>>|<row|<cell|a(v,w)>|<cell|\<assign\>>|<cell|<big|int><rsub|\<Omega\>>\<nabla\>v\<cdot\>\<nabla\>w*d*v<space|1fn><with|mode|text|(energy
      inner product in <with|mode|math|L<rsup|2>(\<Omega\>)>)>>>>>
    </eqnarray*>

    and <with|mode|math|\<\|\|\>v\<\|\|\>=<sqrt|<ip|v|v>>>,
    <with|mode|math|\<interleave\>v\<interleave\>=<sqrt|a(v,v)>>.

    <\equation*>
      H<rsup|1>(\<Omega\>)={v\<in\>L<rsup|2>(\<Omega\>):\<partial\><rsub|1>v,\<partial\><rsub|2>v\<in\>L<rsup|2>(\<Omega\>)}.
    </equation*>
  </definition>

  <\theorem>
    <with|mode|math|H<rsup|1>(\<Omega\>)> is a Hilbert space with
    <with|mode|math|\<\|\|\>u\<\|\|\><rsub|<rsub|1>>=<sqrt|\<\|\|\>u\<\|\|\><rsup|2>+\<interleave\>u\<interleave\><rsup|2>>>.
    <with|mode|math|H<rsup|1><rsub|0>(\<Omega\>)\<assign\>clos(C<rsub|0><rsup|\<infty\>>(\<Omega\>),<l2norm|\<cdot\>><rsub|1>)>
    is a Hilbert space with <with|mode|math|<enorm|\<cdot\>>>.
  </theorem>

  <\definition>
    Let <with|mode|math|\<Omega\>\<subset\>\<bbb-R\><rsup|2>> be a polygonal
    domain, <with|mode|math|h\<gtr\>0> the mesh size parameter. A uniform,
    consistent triangulation <with|mode|math|\<cal-C\><rsub|h>> is a
    decomposition on <with|mode|math|<wide|\<Omega\>|\<bar\>>=<big|cup><rsub|C\<in\>\<cal-C\><rsub|h>><wide|\<Omega\>|\<bar\>><rsub|C>>
    such that

    <big-figure|<postscript|tri-transform.fig|*5/8|*5/8||||>|The triangle
    transformation >

    <\enumerate-alpha>
      <item><with|mode|math|<wide|\<Omega\>|\<bar\>><rsub|C>=T<rsub|C>(<wide|\<Omega\>|^>)>
      with <with|mode|math|<wide|\<Omega\>|^>=conv(<wide|z<rsub|0>|^>,<wide|z<rsub|1>|^>,<wide|z<rsub|2>|^>)>
      with

      <\eqnarray*>
        <tformat|<table|<row|<cell|T<rsub|C>(<wide|z|^>)>|<cell|\<assign\>>|<cell|z<rsub|0>+J<rsub|C><wide|z|^>,>>|<row|<cell|D*T<rsub|C>>|<cell|=>|<cell|J<rsub|C>=<matrix|<tformat|<cwith|1|1|1|1|cell-halign|r>|<table|<row|<cell|z<rsub|1>-z<rsub|0>>|<cell|z<rsub|2>-z<rsub|0>>>>>>,>>>>
      </eqnarray*>

      <item><with|mode|math|det J<rsub|C>\<gtr\>0>

      <\eqnarray*>
        <tformat|<table|<row|<cell|\|J<rsub|C>\|>|<cell|\<leqslant\>>|<cell|C<rsub|\<Omega\>>h,>>|<row|<cell|\|J<rsub|C><rsup|-1>\|>|<cell|\<leqslant\>>|<cell|C<rsub|\<Omega\>>h<rsup|-1>,>>>>
      </eqnarray*>

      with <with|mode|math|\|J\|=sup<rsub|\|<wide|z|^>\|=1>\|J*<wide|z|^>\|>,
      <with|mode|math|\|<wide|z|^>\|=<sqrt|<wide|z|^><rsub|1><rsup|2>+<wide|z<rsub|2>|^><rsup|2>>>,
      <with|mode|math|z<rsub|i>=T<rsub|C><wide|z|^><rsub|i>> (mesh
      uniformness),

      <item><with|mode|math|<wide|\<Omega\>|\<bar\>><rsub|C>\<cap\><wide|\<Omega\>|\<bar\>><rsub|C<rprime|'>>=conv({z<rsub|0>,z<rsub|1>,z<rsub|2>}\<cap\>{z<rsub|0><rprime|'>,z<rsub|1><rprime|'>,z<rsub|2><rprime|'>}>
      (mesh admissibility).
    </enumerate-alpha>
  </definition>

  <\theorem>
    <\enumerate-alpha>
      <item><with|mode|math|X<rsub|h>\<assign\>{v\<in\>C(<wide|\<Omega\>|\<bar\>>):v<rsub|h>\|<rsub|<wide|\<Omega\>|\<bar\>><rsub|C>><with|mode|text|
      linear>}\<subset\>H<rsup|1>(\<Omega\>)>

      <item><with|mode|math|v<rsub|h>\<in\>X<rsub|h>> is uniquely defined by
      the nodal values <with|mode|math|v<rsub|h>(z)> at
      <with|mode|math|z\<in\><wide|\<cal-N\>|\<bar\>><rsub|h>\<assign\><big|cup><rsub|C\<in\>\<cal-C\><rsub|h>><big|cup><rsub|i>T<rsub|C>(<wide|z|^><rsub|i>)>.

      <item><with|mode|math|V<rsub|h>\<assign\>X<rsub|h>\<cap\>H<rsup|1><rsub|0>(\<Omega\>)=span{\<varphi\><rsub|z>\<in\>X<rsub|h>:z\<in\>\<cal-N\><rsub|h>\<assign\><wide|\<cal-N\>|\<bar\>><rsub|h>\<setminus\>\<partial\>\<Omega\>}>
      with

      <\equation*>
        \<varphi\><rsub|z>(y)=<choice|<tformat|<cwith|1|1|1|1|cell-valign|b>|<table|<row|<cell|1>|<cell|y=z>>|<row|<cell|0>|<cell|<with|mode|text|elsewhere>>>>>><space|1fn>(y\<in\>\<cal-N\><rsub|h>).
      </equation*>
    </enumerate-alpha>
  </theorem>

  <\lemma>
    The finite element problem <with|mode|math|u<rsub|h>\<in\>V<rsub|h>:><with|mode|math|a(u<rsub|h>,v<rsub|h>)=(f,v<rsub|h>)>
    for <with|mode|math|v\<in\>\<cal-V\><rsub|h>> has a unique solution.
  </lemma>

  <\proof>
    <with|mode|math|\<b-A\>=(a(\<varphi\><rsub|z>,\<varphi\><rsub|y>))<rsub|z,y\<in\>\<cal-N\><rsub|h>>>,
    <with|mode|math|\<b-f\>=((f,y<rsub|z>))<rsub|z\<in\>\<cal-N\><rsub|h>>>.
    <with|mode|math|\<b-u\>=\<b-A\><rsup|-1>\<b-f\>>, where
    <with|mode|math|u<rsub|h>=<big|sum><rsub|z\<in\>\<cal-N\><rsub|h>>\<b-u\>[z]\<varphi\><rsub|z>>.
    Lax-Milgram ensures the existence of <with|mode|math|A<rsup|-1>>.
  </proof>

  <\example*>
    <with|mode|math|\<Omega\>=(0,1)<rsup|2>>, <with|mode|math|h=1/n>. Then
    <with|mode|math|\<b-A\>> is exactly the same as in Section
    <reference|subsec:fd-model-problem>.

    We obtain <with|mode|math|\|J<rsub|C>\|=h>,
    <with|mode|math|\|J<rsub|C><rsup|-1>\|=h<rsup|-1>><with|mode|math|\<Rightarrow\>><with|mode|math|C<rsub|\<Omega\>>=1>.

    <big-figure|<postscript|uniform-grid-unit-square.fig|4cm|||||>|A uniform
    grid on <with|mode|math|(0,1)<rsup|2>>.>
  </example*>

  We define <with|mode|math|A<rsub|h>:V<rsub|h>\<rightarrow\>V<rsub|h>> to be
  <with|mode|math|(A<rsub|h>v<rsub|h>,w<rsub|h>)=a(v<rsub|h>,w<rsub|h>)> for
  <with|mode|math|v<rsub|h>,w<rsub|h>\<in\>V<rsub|h>>.

  <\equation*>
    \<b-u\><rsup|k+1>=\<b-u\><rsup|k>+\<theta\>(\<b-f\>-\<b-A\>\<b-u\><rsup|k>)
  </equation*>

  corresponds to

  <\equation*>
    u<rsup|k+1><rsub|h>=u<rsup|k><rsub|h>\<theta\><rsub|h>(f<rsub|h>-A<rsub|h>u<rsub|h><rsup|k>)<rsup|h>.
  </equation*>

  The former converges if <with|mode|math|\<rho\>(\<b-I\>-\<theta\>\<b-A\>)=\<\|\|\>\<b-I\>-\<theta\>\<b-A\>\<\|\|\>\<less\>1>.
  (For the equality, remember that <with|mode|math|\<b-A\>> is symmetric.)
  The latter converges if <with|mode|math|\<\|\|\>id-\<theta\><rsub|h>A<rsub|h>\<\|\|\>\<less\>1\<Leftrightarrow\>\<theta\><rsub|h>\<\|\|\>A<rsub|h>\<\|\|\>\<less\>1>,
  where <with|mode|math|\<\|\|\>A<rsub|h>\<\|\|\>\<assign\>sup<rsub|\<\|\|\>v<rsub|h>\<\|\|\>=1>\<\|\|\>A<rsub|h>v<rsub|h>\<\|\|\>>.

  <big-figure|<postscript|red-refinement.fig|*5/8|*5/8||||>|Red refinement>

  Let <with|mode|math|\<cal-C\><rsub|H>> be a uniform, consistent
  triangulation. Then, define <with|mode|math|\<cal-C\><rsub|h>> with
  <with|mode|math|h=1/2H> by uniform refinement, as in

  <\equation*>
    <wide|\<Omega\>|\<bar\>>=<big|cup><rsub|C\<in\>\<cal-C\><rsub|H>><big|cup><rsub|i=0><rsup|3>T<rsub|C>(<wide|\<Omega\>|^><rsub|i>)
  </equation*>

  with

  <\eqnarray*>
    <tformat|<table|<row|<cell|<wide|\<Omega\>|^><rsub|0>>|<cell|=>|<cell|conv{<wide|z|^><rsub|0>,<wide|z|^><rsub|01>,<wide|z|^><rsub|02>},>>|<row|<cell|<wide|\<Omega\>|^><rsub|1>>|<cell|=>|<cell|conv
    {<wide|z|^><rsub|01>,<wide|z|^><rsub|12>,<wide|z|^><rsub|02>},>>|<row|<cell|<wide|\<Omega\>|^><rsub|2>>|<cell|=>|<cell|conv
    {<wide|z|^><rsub|01>,<wide|z|^><rsub|1>,<wide|z|^><rsub|12>},>>|<row|<cell|<wide|\<Omega\>|^><rsub|3>>|<cell|=>|<cell|conv
    {<wide|z|^><rsub|02>,<wide|z|^><rsub|12>,<wide|z|^><rsub|2>},>>>>
  </eqnarray*>

  which implies <with|mode|math|V<rsub|H>\<subset\>V<rsub|h>><emdash>the
  finite element spaces are <em|nested>.

  <\equation*>
    A<rsub|h>u<rsub|h>=f<rsub|h>\<Leftrightarrow\><ip|A<rsub|h>u<rsub|h>|v<rsub|h>>=<ip|f<rsub|h>|v<rsub|h>><space|1fn><with|mode|text|for
    all <with|mode|math|v<rsub|h>\<in\>V<rsub|h>>>.
  </equation*>

  Now define <with|mode|math|c<rsub|H>\<in\>V<rsub|H>> as a form of
  coarse-space residual by requiring

  <\equation*>
    a(c<rsub|H>,v<rsub|H>)=<ip|A<rsub|H>c<rsub|H>|v<rsub|H>><above|=|!><ip|f<rsub|h>-A<rsub|h>u<rsup|k><rsub|h>|v<rsub|H>>=<ip|f<rsub|h>|v<rsub|H>>-a(u<rsup|k>,v<rsub|H>)<space|1fn><with|mode|text|for
    all <with|mode|math|v<rsub|H>\<in\>V<rsub|H>.>>
  </equation*>

  (<with|mode|math|c<rsub|H>> may also be viewed as an element of
  <with|mode|math|V<rsub|h>> through interpolation.)

  <\note>
    <em|(Two-level method)> Let <with|mode|math|u<rsub|h><rsup|0>\<in\>V<rsub|h>>.
    For <with|mode|math|k=1,2,\<ldots\>>

    <\enumerate>
      <item>compute <with|mode|math|u<rsub|h><rsup|k-1/2>=u<rsup|k-1><rsub|h>+\<theta\><rsub|h>(f<rsub|h>-A<rsub|h>u<rsub|h><rsup|k-1>)>,

      <item>compute <with|mode|math|c<rsup|k><rsub|H>\<in\>V<rsub|h>>:
      <with|mode|math|a(c<rsub|H><rsup|k>,v<rsub|H>)=<ip|f<rsub|h>|v<rsub|H>>-a(u<rsub|h><rsup|k-1/2>,v<rsub|H>)>,

      <item>set <with|mode|math|u<rsup|k>=u<rsub|h><rsup|k-1/2>+c<rsub|H><rsup|k>>
    </enumerate>
  </note>

  <\lemma>
    Let

    <\eqnarray*>
      <tformat|<table|<row|<cell|Q<rsub|H>:L<rsup|2>(\<Omega\>)\<rightarrow\>V<rsub|H><with|mode|text|
      be the <with|mode|math|L<rsup|2>>-projection, i.e.
      >>|<cell|>|<cell|<ip|Q<rsub|H>v|w<rsub|H>>=<ip|v|w<rsub|H>>,>>|<row|<cell|P<rsub|H>:H<rsup|1><rsub|0>(\<Omega\>)\<rightarrow\>V<rsub|H><with|mode|text|
      be the Galerkin projection, i.e. >>|<cell|>|<cell|a(P<rsub|H>v,w<rsub|H>)=a(v,w<rsub|H>)>>>>
    </eqnarray*>

    with <with|mode|math|v\<in\>L<rsup|2>(\<Omega\>)\<supset\>V<rsub|h>>,
    <with|mode|math|w<rsub|H>\<in\>V<rsub|H>>,
    <with|mode|math|v\<in\>H<rsup|1><rsub|0>(\<Omega\>)\<supset\>V<rsub|h>>,
    for all <with|mode|math|w<rsub|H>\<in\>V<rsub|H>>. Then, we have for the
    error propagation of the two-level method

    <\equation*>
      e<rsup|k><rsub|h>=<wide*|u<rsup|k><rsub|h>|\<wide-underbrace\>><rsub|=A<rsup|-1><rsub|h>f<rsub|h>>-u<rsub|h>=(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>.
    </equation*>
  </lemma>

  <\proof>
    <\enumerate-alpha>
      <item>

      <\eqnarray*>
        <tformat|<table|<row|<cell|e<rsup|k-1/2>>|<cell|=>|<cell|u<rsup|k-1/2><rsub|h>-u<rsub|h>>>|<row|<cell|>|<cell|=>|<cell|u<rsup|k-1>-u<rsub|h>+\<theta\><rsub|h>(A<rsub|h>u<rsub|h>-A<rsub|h>u<rsub|h><rsup|k-1>)=(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>.>>>>
      </eqnarray*>

      <item>

      <\equation*>
        (A<rsub|H>c<rsub|H><rsup|k>,v<rsub|H>)=a(c<rsub|H><rsup|k>,v<rsub|H>)=a(u<rsub|h>,v<rsub|H>)-a(u<rsub|h><rsup|k-1/2>,v<rsub|H>)=<ip|Q<rsub|H>f<rsub|h>|v<rsub|H>>-<ip|A<rsub|h>u<rsup|k-1/2><rsub|h>|v<rsub|H>>=<ip|Q<rsub|H>f<rsub|h>-Q<rsub|H>A<rsub|h>u<rsup|k-1/2><rsub|h>|v<rsub|H>>
      </equation*>

      Then

      <\equation*>
        A<rsub|H>c<rsub|H><rsup|k>=Q<rsub|H>(<wide*|f<rsub|h>|\<wide-underbrace\>><rsub|=A<rsub|h>u<rsub|h>>-A<rsub|h>u<rsup|k-1/2><rsub|h>)=-Q<rsub|H>A<rsub|h>e<rsup|k-1/2><rsub|h>,
      </equation*>

      so

      <\equation*>
        e<rsup|k><rsub|h>=u<rsup|k><rsub|h>-u<rsub|h>=u<rsup|k-1/2><rsub|h>-u<rsub|h>+c<rsub|H>=e<rsup|k-1/2><rsub|h>+A<rsub|H><rsup|-1>(-Q<rsub|H>A<rsub|h>e<rsup|k-1/2><rsub|h>).
      </equation*>

      <item>For <with|mode|math|w<rsub|H>\<in\>V<rsub|H>>,

      <\equation*>
        <ip|A<rsub|H>P<rsub|H>v<rsub|h>|w<rsub|H>>=a(P<rsub|H>v<rsub|h>,w<rsub|H>)=a(v<rsub|h>,w<rsub|H>)=<ip|A<rsub|h>v<rsub|h>|w<rsub|H<with|color|red|>>>=<ip|Q<rsub|H>A<rsub|h>v<rsub|h>|w<rsub|H>>,
      </equation*>

      so

      <\equation*>
        A<rsub|H>P<rsub|H>=Q<rsub|H>A<rsub|h>\<Rightarrow\>P<rsub|H>=A<rsup|-1><rsub|H>Q<rsub|H>A<rsub|h>.
      </equation*>
    </enumerate-alpha>

    \;
  </proof>

  <\remark>
    We observe that

    <\itemize-dot>
      <item><with|mode|math|Q<rsub|H>> is an orthogonal projection.
      <with|mode|math|Q<rsub|H><rsup|2>=Q<rsub|H>>,
      <with|mode|math|Q<rsub|H>(id-Q<rsub|H>)=0>,
      <with|mode|math|<l2norm|Q<rsub|H>>=1>,
      <with|mode|math|<l2norm|P<rsub|H>>\<gtr\>1>!

      <item><with|mode|math|P<rsub|H>> is an orthogonal projection.
      <with|mode|math|P<rsub|H><rsup|2>=P<rsub|H>>,
      <with|mode|math|P<rsub|H>(id-P<rsub|H>)=0>,
      <with|mode|math|<enorm|P<rsub|H>>=1>.
    </itemize-dot>
  </remark>

  <\theorem>
    <label|the:twolevel>Let

    <\enumerate>
      <item><with|mode|math|<enorm|v<rsub|h>><rsup|2>\<leqslant\>\<theta\><rsub|h><rsup|-1><l2norm|v<rsub|h>><rsup|2>>
      (<with|mode|math|\<Leftrightarrow\>><with|mode|math|\<theta\><rsub|h><l2norm|A<rsub|h>>\<leqslant\>1>),

      <item><with|mode|math|<l2norm|v<rsub|h>-Q<rsub|H>v<rsub|h>><rsup|2>\<leqslant\>C\<theta\><rsub|h><enorm|v<rsub|h>><rsup|2>>
      (<with|mode|math|\<leqslant\>C<l2norm|v<rsub|h>><rsup|2>><with|mode|math|\<Rightarrow\>><with|mode|math|C\<geqslant\>1>).
    </enumerate>

    Then, we have <with|mode|math|<enorm|(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)><rsup|2>\<leqslant\>1-1/C\<less\>1>,
    i.e. <with|mode|math|<enorm|u<rsup|k><rsub|h>-u<rsub|h>>\<leqslant\>(1-1/C)<rsup|k/2><enorm|u<rsub|h><rsup|0>-u<rsub|h>>>.
  </theorem>

  1) and 2) of Theorem <reference|the:twolevel> imply

  <\equation*>
    <l2norm|v<rsub|h>-Q<rsub|H>v<rsub|h>><rsup|2>\<leqslant\>C<l2norm|v<rsub|h>><rsup|2>\<Rightarrow\>C\<geqslant\>1.
  </equation*>

  The following two lemmata provide the (seemingly strong) assumptions of
  this theorem.

  <\lemma>
    <label|lem:2step-energy-estimate>

    <\enumerate-alpha>
      <item><with|mode|math|<l2norm|\<nabla\><wide|v|^>><rsub|<wide|\<Omega\>|^>>\<leqslant\><sqrt|48><l2norm|<wide|v|^>><rsub|<wide|\<Omega\>|^>>>
      for linear <with|mode|math|<wide|v|^>(x)=(1-<wide|x|^><rsub|1>-<wide|x|^><rsub|2>)v<rsub|0>+<wide|x|^><rsub|1>v<rsub|1>+<wide|x|^><rsub|2>v<rsub|2>>.

      <item><with|mode|math|<l2norm|\<nabla\>v<rsub|h>>\<leqslant\>C<rsub|I>h<rsup|-1><l2norm|v<rsub|h>>>
      for <with|mode|math|v<rsub|h>\<in\>X<rsub|h>> with
      <with|mode|math|C<rsub|I>=<sqrt|48>C<rsub|\<Omega\>>>.
    </enumerate-alpha>
  </lemma>

  <with|mode|math|C<rsub|\<Omega\>>> measures the quality of the
  triangulation. Roughly, <with|mode|math|C<rsub|\<Omega\>>\<thickapprox\>1>
  for ``rectangular'' triangles. If using uniform refinement,
  <with|mode|math|C<rsub|\<Omega\>>> is independent of the mesh size.

  <\proof>
    <\enumerate-alpha>
      <item>Using

      <\equation*>
        \<nabla\><wide|v|^>=<matrix|<tformat|<table|<row|<cell|v<rsub|1>-v<rsub|0>>>|<row|<cell|v<rsub|2>-v<rsub|0>>>>>>,
      </equation*>

      and for <with|mode|math|P\<in\>\<bbb-P\><rsub|2>>

      <\equation*>
        <l2norm|P><rsub|<wide|\<Omega\>|^>>=<frac|1|2>\<cdot\><frac|1|3><left|(>P<left|(><frac|1|2>(<wide|z<rsub|0>|^>+<wide|z<rsub|1>|^>)<right|)>+P<left|(><frac|1|2>(<wide|z<rsub|1>|^>+<wide|z<rsub|2>|^>)<right|)>+P<left|(><frac|1|2>(<wide|z<rsub|0>|^>+<wide|z<rsub|2>|^>)<right|)><right|)>
      </equation*>

      and

      <\equation*>
        <frac|1|2>(a+b)<rsup|2>\<leqslant\>a<rsup|2>+b<rsup|2>,
      </equation*>

      we compute

      <\eqnarray*>
        <tformat|<table|<row|<cell|<l2norm|\<nabla\><wide|v|^><rsub|h>><rsub|<wide|\<Omega\>|^>><rsup|2>>|<cell|=>|<cell|<frac|1|2>(v<rsub|1>-v<rsub|0>)<rsup|2>+<frac|1|2>(v<rsub|2>-v<rsub|0>)<rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|2>(v<rsub|1>-v<rsub|2>+v<rsub|2>-v<rsub|0>)<rsup|2>+<frac|1|2>(v<rsub|2>-v<rsub|1>+v<rsub|1>-v<rsub|0>)<rsup|2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|(v<rsub|1>-v<rsub|2>)<rsup|2>+(v<rsub|2>-v<rsub|0>)<rsup|2>+(v<rsub|2>-v<rsub|1>)<rsup|2>+(v<rsub|1>-v<rsub|0>)<rsup|2>>>|<row|<cell|>|<cell|=>|<cell|2(v<rsub|2>-v<rsub|1>)<rsup|2>+(v<rsub|2>-v<rsub|0>)<rsup|2>+(v<rsub|1>-v<rsub|0>)<rsup|2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|2<left|[>(v<rsub|2>-v<rsub|1>)<rsup|2>+(v<rsub|2>-v<rsub|0>)<rsup|2>+(v<rsub|1>-v<rsub|0>)<rsup|2><right|]>>>|<row|<cell|>|<cell|<with|color|red|<above|\<leqslant\>|?>>>|<cell|2\<cdot\>4\<cdot\>6\<cdot\><frac|1|6>\<cdot\><frac|1|4><left|(>(v<rsub|1>+v<rsub|2>)<rsup|2>+(v<rsub|2>+v<rsub|1>)<rsup|2>+(v<rsub|1>+v<rsub|0>)<rsup|2><right|)>>>|<row|<cell|>|<cell|=>|<cell|48<l2norm|<wide|v|^><rsub|h>><rsub|<wide|\<Omega\>|^>><rsup|2>.>>>>
      </eqnarray*>

      <item>On the transformed triangle, with
      <with|mode|math|x=T<rsub|C>(<wide|x|^>)> and considering
      <with|mode|math|D*v(x)=\<nabla\>v(x)<rsup|T>>,

      <\eqnarray*>
        <tformat|<table|<row|<cell|<wide|v|^>(<wide|x|^>)>|<cell|=>|<cell|v<rsub|h>(x)>>|<row|<cell|\<Rightarrow\>D<wide|v|^>(<wide|x|^>)>|<cell|=>|<cell|D(v\<circ\>T<rsub|C>(<wide|x|^>))=D*v(x)D*T<rsub|C>>>|<row|<cell|\<Rightarrow\>\<nabla\><wide|v|^>(x)>|<cell|=>|<cell|D<wide|v|^>(<wide|x|^>)<rsup|T>=(D*T<rsub|C>)<rsup|T>\<nabla\>v(x).>>>>
      </eqnarray*>

      <\eqnarray*>
        <tformat|<table|<row|<cell|<l2norm|\<nabla\>v><rsub|\<Omega\><rsub|C>><rsup|2>=<big|int><rsub|\<Omega\><rsub|C>>\|\<nabla\>v\|<rsup|2>d*x>|<cell|=>|<cell|<big|int><rsub|<wide|\<Omega\>|^>>\|det*J<rsub|C>\|\|\<nabla\>v\<circ\>T\|<rsup|2>d*<wide|x|^>>>|<row|<cell|>|<cell|=>|<cell|<big|int><rsub|<wide|\<Omega\>|^>>\|det*J<rsub|C>\|\|J<rsub|C><rsup|-T>\<nabla\>v\|<rsup|2>d<wide|x|^>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|\|J<rsub|C><rsup|-T>\|<rsup|2>\|det
        J<rsub|C>\|<l2norm|\<nabla\><wide|v|^>><rsup|2><rsub|<wide|\<Omega\>|^>>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|C<rsub|\<Omega\>>h<rsup|-2>48<l2norm|<wide|v|^>><rsup|2><rsub|\<Omega\><rsub|C>>>>>>
      </eqnarray*>

      So,

      <\equation*>
        <enorm|v<rsub|h>><rsup|2>=<big|sum><rsub|C\<in\>\<cal-C\><rsub|h>><l2norm|\<nabla\>v<rsub|h>><rsub|\<Omega\><rsub|C>>\<leqslant\>48C<rsub|\<Omega\>><rsup|2>h<rsup|-2><big|sum><rsub|C\<in\>\<cal-C\><rsub|h>><l2norm|v<rsub|h>><rsup|2><rsub|\<Omega\><rsub|C>>.
      </equation*>
    </enumerate-alpha>

    \;
  </proof>

  <\corollary>
    Let <with|mode|math|A<rsub|h>:V<rsub|h>\<rightarrow\>V<rsub|h>> be
    defined by <with|mode|math|<ip|A<rsub|h>v<rsub|h>|w<rsub|h>>=a(v<rsub|h>,w<rsub|h>)>.
    Then

    <\equation*>
      <l2norm|A<rsub|h>>=sup<rsub|<l2norm|v<rsub|h>>=1><l2norm|A<rsub|h>v<rsub|h>>\<leqslant\>C<rsub|I>h<rsup|-2>.
    </equation*>

    \;
  </corollary>

  <\proof>
    \;

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|A<rsub|h>v<rsub|h>>>|<cell|=>|<cell|sup<rsub|<l2norm|w<rsub|h>>=1><ip|A<rsub|h>v<rsub|h>|w<rsub|h>>=sup<rsub|<l2norm|w<rsub|h>>=1>a(v<rsub|h>,w<rsub|h>)>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|sup<rsub|<l2norm|w<rsub|h>>=1><enorm|v<rsub|h>><enorm|w<rsub|h>>>>|<row|<cell|>|<cell|<below|\<leqslant\>|(<with|mode|text|<reference|lem:2step-energy-estimate>)>>>|<cell|sup<rsub|<l2norm|w<rsub|h>>=1>C<rsub|I>h<rsup|-1><l2norm|v<rsub|h>>C<rsub|I>h<rsup|-1><l2norm|w<rsub|h>>=C<rsub|I>h<rsup|-2><l2norm|v<rsub|h>>.>>>>
    </eqnarray*>
  </proof>

  <\lemma>
    <\enumerate-alpha>
      <item>Let <with|mode|math|<wide|v|^>\<in\>C(<wide|\<Omega\>|^>)> be
      linear on <with|mode|math|<wide|\<Omega\>|^><rsub|i>>, let
      <with|mode|math|<wide|\<Pi\>|^><wide|v|^>> be linear on
      <with|mode|math|<wide|\<Omega\>|^>> with
      <with|mode|math|<wide|\<Pi\>|^><wide|v|^>(<wide|z|^><rsub|i>)=<wide|v|^>(<wide|z|^><rsub|i>)>.
      Then, we have

      <\equation*>
        <l2norm|<wide|v|^>-<wide|\<Pi\>|^><wide|v|^>><rsub|<wide|\<Omega\>|^>>\<leqslant\>3<l2norm|\<nabla\><wide|v|^>><rsub|<wide|\<Omega\>|^>>.
      </equation*>

      <item>

      <\equation*>
        <l2norm|v<rsub|h>-Q<rsub|H>v<rsub|h>>\<leqslant\>2<sqrt|3>C<rsub|\<Omega\>>h<enorm|v<rsub|h>>.
      </equation*>
    </enumerate-alpha>
  </lemma>

  <\proof>
    <\enumerate-alpha>
      <item>Let <with|mode|math|w<rsub|i,j>\<assign\>v<rsub|i,j>-1/2(v<rsub|i>+v<rsub|j>)>.
      Then

      <\eqnarray*>
        <tformat|<table|<row|<cell|<l2norm|<wide|v|^>-<wide|\<Pi\>|^><wide|v|^>><rsub|<wide|\<Omega\>|^>>>|<cell|=>|<cell|<frac|1|4>*<frac|1|6><left|(><left|(><frac|w<rsub|0,2>|2><right|)><rsup|2>+<left|(><frac|w<rsub|0,2>+w<rsub|0,1>|2><right|)><rsup|2>+<left|(><frac|w<rsub|0,1>|2><right|)><rsup|2><right|)>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<frac|1|4><big|sum><rsub|i\<less\>j><left|(><frac|w<rsub|i,j>|2><right|)><rsup|2>\<leqslant\><frac|1|4><big|sum><rsub|i\<less\>j><left|[>(v<rsub|i,j>-v<rsub|i>)<rsup|2>+(v<rsub|i,j>-v<rsub|j>)<rsup|2><right|]>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|3<l2norm|\<nabla\><wide|v|^>><rsup|2>.>>>>
      </eqnarray*>

      <item>

      <\eqnarray*>
        <tformat|<table|<row|<cell|<l2norm|v<rsub|h>-\<Pi\><rsub|H>v<rsub|H>><rsub|\<Omega\><rsub|C>><rsup|2>>|<cell|=>|<cell|<big|int><rsub|\<Omega\><rsub|C>>\|v<rsub|h>-\<Pi\><rsub|H>v<rsub|H>\|<rsup|2>d*x>>|<row|<cell|>|<cell|=>|<cell|<big|int><rsub|<wide|\<Omega\>|^>>\|det*J<rsub|C>\|\|<wide|v|^>-<wide|\<Pi\>|^><wide|v|^>\|d*<wide|x|^>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|3<big|int><rsub|<wide|\<Omega\>|^>>\|det*J<rsub|C>\|<below|\|\<nabla\><wide|v|^>\|<rsup|2><rsub|>|=\|J<rsub|C><rsup|T>\<nabla\>v<rsub|h>\|>d*x>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|3C<rsub|\<Omega\>><rsup|2>h<rsup|2><l2norm|\<nabla\>v<rsub|h>><rsup|2><rsub|\<Omega\><rsub|C>>.>>>>
      </eqnarray*>

      <\eqnarray*>
        <tformat|<table|<row|<cell|<l2norm|v<rsub|h>-Q<rsub|H>v<rsub|h>><rsup|2>>|<cell|=>|<cell|min<rsub|V<rsub|H>\<in\>V<rsub|H>><l2norm|v<rsub|h>-V<rsub|H>><rsup|2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<l2norm|v<rsub|h>-\<Pi\><rsub|H>v<rsub|h>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|C\<in\>\<cal-C\><rsub|H>><l2norm|v<rsub|h>-\<Pi\>V<rsub|h>><rsup|2><rsub|\<Omega\><rsub|C>>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|3C<rsub|\<Omega\>><rsup|2>h<rsup|2><enorm|v<rsub|h>>.>>>>
      </eqnarray*>
    </enumerate-alpha>
  </proof>

  So, <with|mode|math|\<theta\><rsub|h>=h<rsup|2>/C<rsup|2><rsub|I>> and
  <with|mode|math|C*h<rsup|2>/C<rsup|2><rsub|I>=C<rsup|2><rsub|Q>h<rsup|2>><with|mode|math|\<Rightarrow\>><with|mode|math|C=C<rsup|2><rsub|Q>C<rsup|2><rsub|I>=48\<cdot\>12\<cdot\>C<rsup|4><rsub|\<Omega\>>>.
  We obtain <with|mode|math|\<theta\><rsub|h>=h<rsup|2>/<sqrt|48C<rsub|\<Omega\>>>>.
  <with|color|red|Was um alles in der Welt ist <with|mode|math|C<rsub|Q>>?>
  (s. Wieners 12 unten)

  <\proof>
    (of <reference|the:twolevel>) <with|mode|math|e<rsub|h><rsup|k>=(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>>,\ 

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<rho\>>|<cell|=>|<cell|<enorm|(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)>=sup<rsub|<l2norm|v<rsub|h>>><enorm|(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)>>>|<row|<cell|\<Leftrightarrow\><enorm|e<rsup|k><rsub|h>>>|<cell|=>|<cell|<enorm|(id-P<rsub|H>)<wide*|(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>|\<wide-underbrace\>><rsub|e<rsup|k-1/2><rsub|h>>>\<leqslant\>\<rho\><enorm|e<rsup|k-1><rsub|h>>>>>>
    </eqnarray*>

    for all <with|mode|math|e<rsup|k-1><rsub|h>\<in\>V<rsub|h>>.

    <em|1st step.> For all <with|mode|math|v<rsub|H>\<in\>V<rsub|H>>,

    <\eqnarray*>
      <tformat|<table|<row|<cell|a(e<rsup|k><rsub|h>,v<rsub|H>)>|<cell|=>|<cell|a((id-P<rsub|H>)e<rsup|k-1/2><rsub|h>,v<rsub|H>)>>|<row|<cell|>|<cell|=>|<cell|a((id-P<rsub|H>)<rsup|2>e<rsup|k-1/2><rsub|h>,v<rsub|H>)>>|<row|<cell|>|<cell|=>|<cell|a((id-P<rsub|H>)e<rsup|k><rsub|h>,v<rsub|H>)>>|<row|<cell|>|<cell|=>|<cell|a(e<rsup|k><rsub|h>,(id-P<rsub|H>)v<rsub|H>)=a(e<rsup|k><rsub|h>,v<rsub|H>)-a(e<rsup|k><rsub|h>,<wide*|P<rsub|H>v<rsub|H>|\<wide-underbrace\>><rsub|=v<rsub|H>>)=0.>>>>
    </eqnarray*>

    <em|2nd step.>

    <\eqnarray*>
      <tformat|<table|<row|<cell|<enorm|e<rsub|h><rsup|k>><rsup|2>>|<cell|=>|<cell|a(e<rsup|k><rsub|h>,e<rsup|k><rsub|h>)=a(e<rsup|k><rsub|h>,e<rsup|k><rsub|h>-Q<rsub|H>e<rsup|k><rsub|h>)>>|<row|<cell|>|<cell|=>|<cell|<ip|A<rsub|h>e<rsup|k><rsub|h>|e<rsub|h><rsup|k>-Q<rsub|H>e<rsup|k><rsub|h>>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<l2norm|A<rsub|h>e<rsup|k><rsub|h>><l2norm|e<rsub|h><rsup|k>-Q<rsub|H>e<rsup|k><rsub|h>>\<leqslant\><l2norm|A<rsub|h>e<rsub|h><rsup|k>><sqrt|C\<theta\><rsub|h>><enorm|e<rsup|k><rsub|h>>,>>>>
    </eqnarray*>

    using assumption 2.

    <\equation>
      <label|eq:twolevel-step2>\<Rightarrow\><enorm|e<rsup|k><rsub|h>><rsup|2>\<leqslant\><sqrt|C\<theta\><rsub|h>><l2norm|A<rsub|h>e<rsub|h><rsup|k>>\<Rightarrow\><enorm|e<rsub|h><rsup|k>><rsup|2>\<leqslant\>C\<theta\><rsub|h><l2norm|A<rsub|h>e<rsup|k><rsub|h>><rsup|2>
    </equation>

    \;

    <em|3rd step.>

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<Rightarrow\><enorm|(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k><rsub|h>><rsup|2>>|<cell|=>|<cell|<enorm|e<rsub|h><rsup|k>><rsup|2>-2a(e<rsup|k><rsub|h>,\<theta\><rsub|h>A<rsub|h>e<rsup|k><rsub|h>)+<enorm|\<theta\><rsub|h>A<rsub|h>e<rsup|k><rsub|h>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<enorm|e<rsub|h><rsup|k>><rsup|2>-2\<theta\><rsub|h><l2norm|A<rsub|h>e<rsub|h><rsup|k>><rsup|2>+\<theta\><rsub|h><rsup|2><enorm|A<rsub|h>e<rsup|k><rsub|h>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<enorm|e<rsup|k><rsub|h>><rsup|2>-\<theta\><rsub|h><l2norm|A<rsub|h>e<rsub|h><rsup|k>><rsup|2>-\<theta\><rsub|h><left|(><l2norm|A<rsub|h>e<rsup|k><rsub|h>><rsup|2>-\<theta\><rsub|h><enorm|A<rsub|h>e<rsup|h><rsub|k>><rsup|2><right|)>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<enorm|e<rsup|k><rsub|h>><rsup|2>-\<theta\><rsub|h><l2norm|A<rsub|h>e<rsup|k><rsub|h>><rsup|2>>>|<row|<cell|>|<cell|<below|\<leqslant\>|(<reference|eq:twolevel-step2>)>>|<cell|<enorm|e<rsup|k><rsub|h>><rsup|2>-<frac|1|C><enorm|e<rsup|k><rsub|h>><rsup|2>=<left|(>1-<frac|1|C><right|)><enorm|e<rsup|k><rsub|h>><rsup|2>,>>>>
    </eqnarray*>

    where, considering (<reference|eq:twolevel-step2>), we used

    <\equation*>
      -\<theta\><rsub|h><l2norm|A<rsub|h>e<rsup|k><rsub|h>><rsup|2>\<leqslant\>-<frac|1|C><enorm|e<rsup|k><rsub|h>>.
    </equation*>

    <em|4th step.>

    <\eqnarray*>
      <tformat|<table|<row|<cell|<enorm|e<rsup|k><rsub|h>><rsup|2>>|<cell|=>|<cell|<enorm|(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|a((id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>,(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>)>>|<row|<cell|>|<cell|=>|<cell|a((id-\<theta\><rsub|h>A<rsub|h>)(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>,e<rsup|k-1><rsub|h>)>>|<row|<cell|>|<cell|=>|<cell|<ip|(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k><rsub|h>|e<rsup|k><rsub|h>>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<enorm|(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k><rsub|h>><enorm|e<rsup|k-1><rsub|h>>\<leqslant\><rsub|3rd
      step><sqrt|1-1/c><enorm|e<rsup|k><rsub|h>><enorm|e<rsup|k-1><rsub|h>>.(id-\<theta\><rsub|h>A<rsub|h>)>>>>
    </eqnarray*>

    \;
  </proof>

  <subsection|Implementation of the Two-Level Method>

  <\eqnarray*>
    <tformat|<table|<row|<cell|A<rsub|H>>|<cell|:>|<cell|V<rsub|H>\<rightarrow\>V<rsub|H>>>|<row|<cell|>|<cell|>|<cell|(A<rsub|H>v<rsub|H>,w<rsub|H>)=a(v<rsub|H>,w<rsub|H>)>>|<row|<cell|\<b-A\><rsub|H>>|<cell|\<assign\>>|<cell|(a(\<varphi\><rsub|z><rsup|H>,\<varphi\><rsub|y><rsup|H>))<rsub|z,y\<in\>\<cal-N\><rsub|H>>\<in\>\<bbb-R\><rsup|N<rsub|H>\<times\>N<rsub|H>><with|mode|text|
    with >N<rsub|H>:=\|\<cal-N\><rsub|H>\|>>|<row|<cell|>|<cell|>|<cell|\<b-v\><rsub|H><rsup|T>\<b-A\><rsub|H>\<b-w\><rsub|H>=a(v<rsub|H>,w<rsub|H>)>>|<row|<cell|<with|mode|text|for
    >\<b-v\><rsub|H>>|<cell|=>|<cell|<big|sum>\<b-v\><rsub|H>[z]\<varphi\><rsub|z><rsup|H>\<Leftrightarrow\>\<b-v\><rsub|H>=(v<rsub|H>(z))<rsub|z\<in\>\<cal-N\><rsub|H>>>>|<row|<cell|\<b-A\><rsub|h>>|<cell|\<in\>>|<cell|\<bbb-R\><rsup|N<rsub|h>\<times\>N<rsub|h>>>>|<row|<cell|\<b-v\><rsub|H>>|<cell|\<in\>>|<cell|V<rsub|H>\<subset\>V<rsub|h>>>|<row|<cell|\<b-v\><rsub|h>>|<cell|\<assign\>>|<cell|(v<rsub|H>(z))<rsub|z\<in\>\<cal-N\><rsub|h>>>>|<row|<cell|>|<cell|>|<cell|\<Rightarrow\>\<b-v\><rsub|n>(z)=<frac|1|2><left|(>\<b-v\><rsub|h>(x)+\<b-v\><rsub|h>(y)<right|)>>>|<row|<cell|\<b-v\><rsub|h>>|<cell|=>|<cell|\<b-I\><rsub|h><with|mode|text|
    with >\<b-I\><rsub|h>\<in\>\<bbb-R\><rsup|N<rsub|h>\<times\>N<rsub|H>>>>>>
  </eqnarray*>

  <big-figure|<postscript|2step-refinement.fig|7cm|||||>|Interpolation to a
  finer discretization.>

  <\eqnarray*>
    <tformat|<table|<row|<cell|A<rsub|h>u<rsub|h>=f<rsub|h>>|<cell|\<Leftrightarrow\>>|<cell|<ip|A<rsub|h>u<rsub|h>|\<varphi\><rsub|z><rsup|h>>=<ip|f<rsub|h>|\<varphi\><rsub|z><rsup|h>><space|1fn>\<forall\>z\<in\>\<cal-N\><rsub|h>>>|<row|<cell|>|<cell|\<Leftrightarrow\>>|<cell|\<b-A\><rsub|h>\<b-u\><rsub|h>=\<b-f\><rsub|h><with|mode|text|
    with >\<b-f\><rsub|h>=(<ip|f<rsub|h>|\<varphi\><rsub|h><rsup|z>>)<rsub|z\<in\>\<cal-N\><rsub|h>>>>>>
  </eqnarray*>

  Within finite element context,

  <\eqnarray*>
    <tformat|<table|<row|<cell|u<rsub|h><rsup|k-1/2>>|<cell|=>|<cell|u<rsub|h><rsup|k-1>+\<theta\><rsub|h>(f<rsub|h>-A<rsub|h>u<rsub|h><rsup|k-1>)>>|<row|<cell|<with|mode|text|solve
    >a(c<rsub|H>,\<varphi\><rsub|z><rsup|H>)>|<cell|=>|<cell|<ip|f<rsub|h>|\<varphi\><rsub|z><rsup|H>>-a(u<rsup|k-1/2><rsub|h>,\<varphi\><rsub|z><rsup|H>)>>|<row|<cell|<ip|A<rsub|H>c<rsub|H>|\<varphi\><rsub|z><rsup|H>>>|<cell|=>|<cell|<ip|f<rsub|h>-A<rsub|h>u<rsup|k-1/2><rsub|h>|\<varphi\><rsub|z><rsup|H>>.>>>>
  </eqnarray*>

  Within vector context,

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-u\><rsub|h><rsup|k-1/2>>|<cell|=>|<cell|\<theta\><rsub|h>(\<b-f\><rsub|h>-\<b-A\><rsub|h>\<b-u\><rsub|h><rsup|k-1>)>>|<row|<cell|(\<b-A\><rsub|H>\<b-c\><rsub|H>)<rsup|T>\<b-v\><rsub|H>>|<cell|=>|<cell|(\<b-f\><rsub|h>-\<b-A\><rsub|h>\<b-u\><rsub|h><rsup|k-1/2>)(\<b-I\><rsub|h>\<b-v\><rsub|H>)>>|<row|<cell|\<Rightarrow\>\<b-v\><rsup|T><rsub|H>\<b-A\><rsub|H>\<b-c\><rsub|H>>|<cell|=>|<cell|\<b-v\><rsub|H><rsup|T>\<b-I\><rsub|h><rsup|T>(\<b-f\><rsub|h>-\<b-A\><rsub|h>\<b-u\><rsub|h><rsup|k-1/2>)>>|<row|<cell|\<Rightarrow\>\<b-A\><rsub|H>\<b-c\><rsub|H>>|<cell|=>|<cell|\<b-I\><rsub|h><rsup|T>(\<b-f\><rsub|h>-\<b-A\><rsub|h>\<b-u\><rsub|h><rsup|k-1/2>).>>>>
  </eqnarray*>

  In particular,

  <\eqnarray*>
    <tformat|<table|<row|<cell|<ip|A<rsub|h>v<rsub|H>|w<rsub|H>>>|<cell|=>|<cell|<ip|A<rsub|H>v<rsub|H>|w<rsub|H>>>>|<row|<cell|\<Leftrightarrow\>(\<b-A\><rsub|h>\<b-I\><rsub|h>\<b-v\><rsub|H>)<rsup|T>(\<b-I\>\<b-w\><rsub|H>)>|<cell|=>|<cell|(\<b-A\><rsub|H>\<b-v\><rsub|H>)<rsup|T>\<b-w\><rsub|H>>>|<row|<cell|(\<b-v\><rsub|H>)<rsup|T>\<b-I\><rsub|h><rsup|T>>|<cell|=>|<cell|\<b-A\><rsub|h>\<b-I\><rsub|h>\<b-w\><rsub|H>>>|<row|<cell|\<Rightarrow\>\<b-A\><rsub|H>>|<cell|=>|<cell|\<b-I\><rsub|h><rsup|T>\<b-A\><rsub|h>\<b-I\><rsub|h><with|mode|text|
    (called the <em|Galerkin product>)>.>>>>
  </eqnarray*>

  <section|Classical two-level analysis>

  Let <with|mode|math|H<rsub|0>\<supset\>H<rsub|1>\<supset\>H<rsub|2>> be a
  sequence of Hilbert spaces with norms <with|mode|math|<l2norm|\<cdot\>><rsub|m>>.
  Consider a bilinear form

  <\eqnarray*>
    <tformat|<table|<row|<cell|a>|<cell|:>|<cell|H<rsub|1>\<times\>H<rsub|1>\<rightarrow\>\<bbb-R\>>>|<row|<cell|<with|mode|text|which
    is continuous >a(v,w)>|<cell|\<leqslant\>>|<cell|C<rsub|1><l2norm|v><rsub|1><l2norm|w><rsub|1>>>|<row|<cell|<with|mode|text|and
    elliptic >a(v,v)>|<cell|\<geqslant\>>|<cell|C<rsub|0><l2norm|v><rsub|1><rsup|2>.>>>>
  </eqnarray*>

  If <with|mode|math|a(\<cdot\>,\<cdot\>)> is symmetric, set
  <with|mode|math|<enorm|\<cdot\>>=<sqrt|a(v,v)>>. We assume full regularity:
  <with|mode|math|f\<in\>H<rsub|0>> and <with|mode|math|u\<in\>H<rsub|1>>
  with <with|mode|math|a(v,u)=<ip|f|v><rsub|0>> for all
  <with|mode|math|v\<in\>H<rsub|1>><with|mode|math|\<Rightarrow\>><with|mode|math|u\<in\>H<rsub|2>>
  and <with|mode|math|<l2norm|u><rsub|2>\<leqslant\>C<rsub|2><l2norm|f><rsub|0>>.
  (<with|color|red|Was bedeutet der Spaÿ?>)

  Because of this assumption of full regularity, multigrid may not be a great
  idea for non-smooth problems, e.g. hyperbolic equations.

  <\example>
    <with|mode|math|H<rsub|0>\<assign\>L<rsup|2>(0,1)>,
    <with|mode|math|H<rsub|1>\<assign\>{v\<in\>H<rsub|0>:v<rprime|'>\<in\>L<rsup|2>(0,1),v(0)=v(1)=0}>,
    <with|mode|math|H<rsub|2>\<assign\>{v\<in\>H<rsub|1>:v<rprime|''>\<in\>L<rsup|2>(0,1)}>.

    <\equation*>
      a(v,w)\<assign\><big|int><rsub|0><rsup|1>v<rprime|'>w<rprime|'>d*x=-<big|int><rsub|0><rsup|1>v<rprime|''>w*d*x
    </equation*>

    For norms, use

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|v><rsub|0>>|<cell|\<assign\>>|<cell|<sqrt|<big|int>v<rsup|2>d*x>,>>|<row|<cell|<l2norm|v><rsub|1>>|<cell|=>|<cell|<l2norm|v<rprime|'>><rsub|0>,>>|<row|<cell|<l2norm|v><rsub|2>>|<cell|=>|<cell|<l2norm|v<rprime|''>><rsub|0>.>>>>
    </eqnarray*>

    Poincaré-Friedrichs says that <with|mode|math|<l2norm|\<cdot\>><rsub|1>>
    is a norm.

    <\equation*>
      a(u,v)=<ip|f|v>
    </equation*>

    Then <with|mode|math|<l2norm|u<rprime|''>><rsub|0>=<l2norm|u><rsub|2>\<leqslant\><l2norm|f><rsub|0>>,
    and, really, <with|mode|math|C<rsub|0>=C<rsub|1>=C<rsub|2>=1>.
  </example>

  <\example>
    <with|mode|math|\<Omega\>\<in\>\<bbb-R\><rsup|d>> is our domain.

    <\eqnarray*>
      <tformat|<table|<row|<cell|H<rsub|0>>|<cell|\<assign\>>|<cell|L<rsup|2>(\<Omega\>)>>|<row|<cell|H<rsub|1>>|<cell|\<assign\>>|<cell|{v\<in\>H<rsub|0>:\<nabla\>v\<in\>L<rsup|2>(\<Omega\>),v\|<rsub|\<Gamma\>>=0}<with|mode|text|
      with >\<Gamma\>\<subset\>\<partial\>\<Omega\>,
      vol<rsub|d-1>(\<Gamma\>)\<gtr\>0>>|<row|<cell|H<rsub|2>>|<cell|\<assign\>>|<cell|{v\<in\>H<rsub|1>:D<rsub|2>v\<in\>L<rsup|2>(\<Omega\>)<rsup|d\<times\>d>}>>>>
    </eqnarray*>

    Then

    <\equation*>
      a(v,w)=<big|int><rsub|\<Omega\>><left|[>(K(x)\<nabla\>v(x))\<cdot\>\<nabla\>w(x)+c(x)\<cdot\>\<nabla\>v(x)w(x)+r(x)v(x)w(x)<right|]>dx
    </equation*>

    where

    <\eqnarray*>
      <tformat|<table|<row|<cell|K>|<cell|\<in\>>|<cell|L<rsup|\<infty\>>(\<Omega\>)<rsup|d\<times\>d><with|mode|text|
      with >\<xi\><rsup|T>K(x)\<xi\>\<geqslant\>k<rsub|0>\|\<xi\>\|<rsup|2>,
      k<rsub|0>\<gtr\>0<space|1fn>(\<xi\>\<in\>\<bbb-R\><rsup|d>),>>|<row|<cell|c>|<cell|\<in\>>|<cell|L<rsup|\<infty\>>(\<Omega\>)<rsup|d><with|mode|text|
      with >\<nabla\>\<cdot\>c=div*c\<in\>L<rsup|\<infty\>>(\<Omega\>)<with|mode|text|
      almost everywhere in <with|mode|math|\<Omega\>>>,>>|<row|<cell|r>|<cell|\<in\>>|<cell|L<rsup|\<infty\>>(\<Omega\>)<with|mode|text|
      with >r(x)-1/2div*c(x)\<geqslant\>r<rsub|0>\<geqslant\>0<with|mode|text|
      almost everywhere in <with|mode|math|\<Omega\>>>.>>>>
    </eqnarray*>

    Our previous special case was <with|mode|math|K=id>,
    <with|mode|math|r=0>, <with|mode|math|c=0>, with
    <with|mode|math|<enorm|v>=<l2norm|\<nabla\>v>> a norm.

    The Poincaré inequality gives

    <\equation*>
      <l2norm|v><rsub|0>\<leqslant\>C<rsub|p><l2norm|\<nabla\>v><rsub|0><with|mode|text|
      with >C<rsub|p><with|mode|text| dependent on >\<Gamma\>.
    </equation*>

    We obtain

    <\equation*>
      <l2norm|\<nabla\>v>=<enorm|v>\<geqslant\>(1+C<rsub|p><rsup|2>)<rsup|-1/2><l2norm|v><rsub|1>,
    </equation*>

    using

    <\equation*>
      <l2norm|v><rsub|1>=<sqrt|<l2norm|v><rsub|0><rsup|2>+<l2norm|\<nabla\>v><rsup|2><rsub|0>>.
    </equation*>

    In our special case of <with|mode|math|K=id>, <with|mode|math|r=0>,
    <with|mode|math|c=0> and <with|mode|math|\<Omega\>> convex,

    <\equation*>
      <l2norm|D<rsup|2>u><rsub|0>\<leqslant\><wide*|<l2norm|\<Delta\>u><rsub|0>|\<wide-underbrace\>><rsub|=tr
      D<rsup|2>u>\<Rightarrow\><l2norm|><rsub|2>\<leqslant\>C<rsub|2><l2norm|f><rsub|0>.
    </equation*>
  </example>

  <\example>
    <with|mode|math|\<Omega\>\<in\>\<bbb-R\><rsup|d>> is our domain.

    <\eqnarray*>
      <tformat|<table|<row|<cell|H<rsub|0>>|<cell|\<assign\>>|<cell|L<rsup|2>(\<Omega\>)<rsup|d>>>|<row|<cell|H<rsub|1>>|<cell|\<assign\>>|<cell|{v\<in\>H<rsub|0>:D*v\<in\>L<rsup|2>(\<Omega\>),v\|<rsub|\<Gamma\>>=0}<with|mode|text|
      with >\<Gamma\>\<subset\>\<partial\>\<Omega\>,
      vol<rsub|d-1>(\<Gamma\>)\<gtr\>0>>|<row|<cell|H<rsub|2>>|<cell|\<assign\>>|<cell|{v\<in\>H<rsub|1>:D<rsub|2>v<rsub|j>\<in\>L<rsup|2>(\<Omega\>)<rsup|d\<times\>d>}>>>>
    </eqnarray*>

    <with|mode|math|G,F\<in\>\<bbb-R\><rsup|d\<times\>d>>,
    <with|mode|math|sym(F)=1/2(F+F<rsup|T>)>,
    <with|mode|math|F:G=<big|sum>F<rsub|i,j>G<rsub|i,j>>

    <\equation*>
      a(v,w)=<big|int><rsub|\<Omega\>>(2\<mu\>*sym*\<nabla\>v:sym*\<nabla\>w+\<lambda\>div*v*div*w)d*x
    </equation*>

    Korn's inequality gives

    <\equation*>
      a(v,v)\<geqslant\>2\<mu\><l2norm|sym*\<nabla\>v>\<geqslant\>2\<mu\>C<rsub|K><l2norm|\<nabla\>v>
    </equation*>

    with <with|mode|math|\<lambda\>\<geqslant\>0>,
    <with|mode|math|\<mu\>\<gtr\>0>.
  </example>

  <\lemma>
    <\equation*>
      <l2norm|(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)<rsup|m>><rsub|0>\<leqslant\>C,
    </equation*>

    if

    <\enumerate-alpha>
      <item><with|mode|math|<l2norm|(id-P<rsub|H>)v<rsub|h>><rsub|0>\<leqslant\><wide|C|~><rsub|h><l2norm|A<rsub|h>v<rsub|h>><rsub|0>>
      for <with|mode|math|v<rsub|h>\<in\>V<rsub|h>>, <em|(approximation
      property)>

      <item><with|mode|math|<wide|C|~><rsub|h><l2norm|A<rsub|h>(id-\<theta\><rsub|h>A<rsub|h>)e<rsub|h>>\<leqslant\>C<l2norm|e<rsub|h>>>
      for <with|mode|math|e<rsub|h>\<in\>V<rsub|h>>. <em|(smoothing
      property)>
    </enumerate-alpha>
  </lemma>

  <\proof>
    <with|mode|math|v<rsub|h>=(id-\<theta\><rsub|h>A<rsub|h>)<rsup|m>e<rsub|h>>
    and <with|mode|math|<l2norm|(id-P<rsub|H>)(id-Q<rsub|H>A<rsub|h>)<rsup|m>>=sup<rsub|<l2norm|e<rsub|h>><rsub|0>=1><l2norm|(id-P<rsub|H>)(id-Q<rsub|H>A<rsub|h>)<rsup|m>e<rsub|h>><rsub|0>>.
  </proof>

  <\theorem>
    <with|mode|math|V<rsub|H>\<subset\>V<rsub|h>\<subset\>H<rsub|1>> be
    conforming finite element spaces with mesh size <with|mode|math|h>,
    <with|mode|math|H=2h> and

    <\enumerate-alpha>
      <item><with|mode|math|inf<rsub|v<rsub|H>\<in\>V<rsub|H>><l2norm|v-v<rsub|H>>\<leqslant\>C<rsub|Q>H<l2norm|v><rsub|2>>
      for <with|mode|math|v\<in\>H<rsub|2>>, (<with|color|red|sup??>)

      <item><with|mode|math|<l2norm|v<rsub|h>><rsub|1>\<leqslant\>C<rsub|I>h<rsup|-1><l2norm|v<rsub|h>><rsub|0>>
      for <with|mode|math|v<rsub|h>\<in\>V<rsub|h>>.
    </enumerate-alpha>

    Then, the <em|approximation property>

    <\equation*>
      <l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>>\<leqslant\>C<rsub|A>\<lambda\><rsub|h><rsup|-1><l2norm|A<rsub|h>v<rsub|h>><rsub|0>
    </equation*>

    for all <with|mode|math|v<rsub|h>\<in\>V<rsub|h>> holds with

    <\equation*>
      \<lambda\><rsub|n>=<l2norm|A<rsub|h>><rsub|0>=sup<rsub|<l2norm|v<rsub|h>>=1><l2norm|A*v<rsub|h>><rsub|0>.
    </equation*>
  </theorem>
</body>

<\initial>
  <\collection>
    <associate|language|american>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|<uninit>|2>>
    <associate|auto-10|<tuple|2.3|5>>
    <associate|auto-11|<tuple|2.2|?>>
    <associate|auto-12|<tuple|2.4|?>>
    <associate|auto-13|<tuple|3|?>>
    <associate|auto-2|<tuple|1|2>>
    <associate|auto-3|<tuple|1.1|2>>
    <associate|auto-4|<tuple|1.2|3>>
    <associate|auto-5|<tuple|1.3|3>>
    <associate|auto-6|<tuple|2|4>>
    <associate|auto-7|<tuple|2.1|4>>
    <associate|auto-8|<tuple|2.1|4>>
    <associate|auto-9|<tuple|2.2|5>>
    <associate|eq:twolevel-step2|<tuple|2.1|?>>
    <associate|lem:2step-energy-estimate|<tuple|2.10|?>>
    <associate|subsec:fd-model-probelm|<tuple|1.1|?>>
    <associate|subsec:fd-model-problem|<tuple|1.1|2>>
    <associate|the:twolevel|<tuple|2.9|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|The triangle transformation |<pageref|auto-8>>

      <tuple|normal|A uniform grid on <with|mode|<quote|math>|(0,1)<rsup|2>>.|<pageref|auto-9>>

      <tuple|normal|Red refinement|<pageref|auto-10>>

      <tuple|normal||<pageref|auto-12>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Table
      of contents> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      1<space|2spc>Introduction <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>

      <with|par-left|<quote|1.5fn>|1.1<space|2spc>The model problem
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1.5fn>|1.2<space|2spc>Iterative Methods for
      Linear Systems <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|1.5fn>|1.3<space|2spc>The Multigrid Idea
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      2<space|2spc>A two-level method <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>

      <with|par-left|<quote|1.5fn>|2.1<space|2spc>The model problem
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <with|par-left|<quote|1.5fn>|2.2<space|2spc>Implementation of the
      two-level method <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      3<space|2spc>Classical two-level analysis
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>
    </associate>
  </collection>
</auxiliary>