<TeXmacs|1.0.5>

<project|multigrid.tm>

<style|<tuple|article|number-long-article|skript>>

<\body>
  <assign|section-nr|2><section|Classical Two-Level Analysis>

  Let <with|mode|math|H<rsub|0>\<supset\>H<rsub|1>\<supset\>H<rsub|2>> be a
  sequence of Hilbert spaces with norms <with|mode|math|<l2norm|\<cdot\>><rsub|m>>.
  Consider a bilinear form <with|mode|math|a:H<rsub|1>\<times\>H<rsub|1>\<rightarrow\>\<bbb-R\>>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<with|mode|text|which is continuous
    >a(v,w)>|<cell|\<leqslant\>>|<cell|C<rsub|1><l2norm|v><rsub|1><l2norm|w><rsub|1>>>|<row|<cell|<with|mode|text|and
    elliptic >a(v,v)>|<cell|\<geqslant\>>|<cell|C<rsub|0><l2norm|v><rsub|1><rsup|2>.>>>>
  </eqnarray*>

  If <with|mode|math|a(\<cdot\>,\<cdot\>)> is symmetric, set
  <with|mode|math|<enorm|v>\<assign\><sqrt|a(v,v)>>. We assume full
  regularity: If we have <with|mode|math|f\<in\>H<rsub|0>> and
  <with|mode|math|u\<in\>H<rsub|1>> with <with|mode|math|a(u,v)=<ip|f|v><rsub|0>>
  for all <with|mode|math|v\<in\>H<rsub|1>>, then
  <with|mode|math|u\<in\>H<rsub|2>> and <with|mode|math|<l2norm|u><rsub|2>\<leqslant\>C<rsub|2><l2norm|f><rsub|0>>.

  Because of this assumption of full regularity, multigrid may not be a great
  idea for non-smooth problems, e.g. hyperbolic equations.

  <\example>
    <with|mode|math|H<rsub|0>\<assign\>L<rsup|2>(0,1)>,
    <with|mode|math|H<rsub|1>\<assign\>{v\<in\>H<rsub|0>,v<rprime|'>\<in\>L<rsup|2>(0,1),v(0)=v(1)=0}>,
    <with|mode|math|H<rsub|2>\<assign\>{v\<in\>H<rsub|1>,v<rprime|''>\<in\>L<rsup|2>(0,1)}>.

    <\equation*>
      a(v,w)\<assign\><big|int><rsub|0><rsup|1>v<rprime|'>w<rprime|'>d
      x=-<big|int><rsub|0><rsup|1>v<rprime|''>w*d x.
    </equation*>

    For norms, use

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|v><rsub|0>>|<cell|\<assign\>>|<cell|<sqrt|<big|int><rsub|0><rsup|1>v<rsup|2>d
      x>,>>|<row|<cell|<l2norm|v><rsub|1>>|<cell|\<assign\>>|<cell|<l2norm|v<rprime|'>><rsub|0>,>>|<row|<cell|<l2norm|v><rsub|2>>|<cell|\<assign\>>|<cell|<l2norm|v<rprime|''>><rsub|0>.>>>>
    </eqnarray*>

    Poincaré-Friedrichs says that <with|mode|math|<l2norm|\<cdot\>><rsub|1>>
    is a norm. Consider

    <\equation*>
      a(u,v)=<ip|f|v>.
    </equation*>

    Then <with|mode|math|<l2norm|u<rprime|''>><rsub|0>=<l2norm|u><rsub|2>\<leqslant\><l2norm|f><rsub|0>>,
    and, really, <with|mode|math|C<rsub|0>=C<rsub|1>=C<rsub|2>=1>.
  </example>

  <\example>
    <with|mode|math|\<Omega\>\<in\>\<bbb-R\><rsup|d>> is our domain.

    <\eqnarray*>
      <tformat|<table|<row|<cell|H<rsub|0>>|<cell|\<assign\>>|<cell|L<rsup|2>(\<Omega\>),>>|<row|<cell|H<rsub|1>>|<cell|\<assign\>>|<cell|{v\<in\>H<rsub|0>,\<nabla\>v\<in\>L<rsup|2>(\<Omega\>),v\|<rsub|\<Gamma\>>=0}<with|mode|text|
      with >\<Gamma\>\<subset\>\<partial\>\<Omega\>,
      vol<rsub|d-1>(\<Gamma\>)\<gtr\>0,>>|<row|<cell|H<rsub|2>>|<cell|\<assign\>>|<cell|{v\<in\>H<rsub|1>,D<rsub|2>v\<in\>L<rsup|2>(\<Omega\>)<rsup|d\<times\>d>}.>>>>
    </eqnarray*>

    <\equation*>
      a(v,w)\<assign\><big|int><rsub|\<Omega\>>(K(x)\<nabla\>v(x))\<cdot\>\<nabla\>w(x)+c(x)\<cdot\>\<nabla\>v(x)w(x)+r(x)v(x)w(x)d
      x,
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

    The Poincaré-Friedrichs inequality gives

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
      <l2norm|v><rsub|1>\<assign\><sqrt|<l2norm|v><rsub|0><rsup|2>+<l2norm|\<nabla\>v><rsup|2><rsub|0>>.
    </equation*>

    In our special case of <with|mode|math|K=id>, <with|mode|math|r=0>,
    <with|mode|math|c=0> and <with|mode|math|\<Omega\>> convex,

    <\equation*>
      <l2norm|D<rsup|2>u><rsub|0>\<leqslant\><wide*|<l2norm|\<Delta\>u><rsub|0>|\<wide-underbrace\>><rsub|=tr
      D<rsup|2>u>\<Rightarrow\><l2norm|u><rsub|2>\<leqslant\>C<rsub|2><l2norm|f><rsub|0>.
    </equation*>
  </example>

  <\example>
    <with|mode|math|\<Omega\>\<in\>\<bbb-R\><rsup|d>> is our domain.

    <\eqnarray*>
      <tformat|<table|<row|<cell|H<rsub|0>>|<cell|\<assign\>>|<cell|L<rsup|2>(\<Omega\>)<rsup|d>,>>|<row|<cell|H<rsub|1>>|<cell|\<assign\>>|<cell|{v\<in\>H<rsub|0>,D*v\<in\>L<rsup|2>(\<Omega\>),v\|<rsub|\<Gamma\>>=0}<with|mode|text|
      with >\<Gamma\>\<subset\>\<partial\>\<Omega\>,
      vol<rsub|d-1>(\<Gamma\>)\<gtr\>0,>>|<row|<cell|H<rsub|2>>|<cell|\<assign\>>|<cell|{v\<in\>H<rsub|1>,D<rsub|2>v<rsub|j>\<in\>L<rsup|2>(\<Omega\>)<rsup|d\<times\>d>}.>>>>
    </eqnarray*>

    <with|mode|math|G,F\<in\>\<bbb-R\><rsup|d\<times\>d>>,
    <with|mode|math|sym(F)=1/2(F+F<rsup|T>)>,
    <with|mode|math|F:G=<big|sum>F<rsub|i,j>G<rsub|i,j>>

    <\equation*>
      a(v,w)\<assign\><big|int><rsub|\<Omega\>>(2\<mu\>*sym*\<nabla\>v:sym*\<nabla\>w+\<lambda\>div*v*div*w)d
      x.
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

  <\remark*>
    \;

    <\eqnarray*>
      <tformat|<table|<row|<cell|a(u,v)>|<cell|=>|<cell|<ip|f|v><space|1em>u\<in\>H<rsub|2>,v\<in\>H<rsup|1>>>|<row|<cell|a(u<rsub|h>,v<rsub|h>)>|<cell|=>|<cell|<ip|f|f<rsub|h>><space|1em>v<rsub|h>\<in\>V<rsub|h>>>>>
    </eqnarray*>

    By standard finite element analysis, we find

    <\equation*>
      <l2norm|u-u<rsub|h>><rsub|1>\<leqslant\>C<rsub|3>inf<rsub|v<rsub|h>\<in\>V<rsub|h>><l2norm|u-v<rsub|h>>\<leqslant\>C<rsub|4>h<l2norm|u><rsub|2>\<leqslant\>C<rsub|5><wide*|h|\<wide-underbrace\>><rsub|\<lambda\><rsub|h>=h<rsup|-2>><wide*|<l2norm|f><rsub|0>|\<wide-underbrace\>><rsub|=A<rsub|h>u<rsub|h>\<Rightarrow\><l2norm|u-u<rsub|h>><rsub|0>\<leqslant\>C<rsub|6>h<rsup|2><l2norm|f><rsub|0>>+<with|mode|text|duality>.
    </equation*>
  </remark*>

  <\proof>
    <em|1st step.> Consider a dual solution <with|mode|math|z\<in\>H<rsub|1>>
    of

    <\equation*>
      a(v,z)=<ip|v<rsub|h>-P<rsub|H>v<rsub|h>|v><rsub|0><space|1em>\<forall\>v\<in\>H<rsub|1>.
    </equation*>

    By regularity, we get <with|mode|math|z\<in\>H<rsub|2>\<Rightarrow\><l2norm|z><rsub|2>\<leqslant\>C<rsub|2><l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|0>>.

    Since <with|mode|math|P<rsub|H>(id-P<rsub|H>)=0>, we have
    <with|mode|math|a(v<rsub|h>-P<rsub|H>v<rsub|h>,w<rsub|H>)=0> and

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|0><rsup|2>>|<cell|=>|<cell|a(v<rsub|h>-P<rsub|H>v<rsub|h>,z)=a(v<rsub|h>-P<rsub|H>v<rsub|h>,z-w<rsub|H>)\<leqslant\>C<rsub|1><l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><l2norm|z-w<rsub|H>>>>|<row|<cell|\<Rightarrow\><l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|0><rsup|2>>|<cell|\<leqslant\>>|<cell|C<rsub|1>C<rsub|Q>H<l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|1><l2norm|z><rsub|2>\<leqslant\>C<rsub|1>C<rsub|2>C<rsub|Q><wide*|H|\<wide-underbrace\>><rsub|=2h><l2norm|v<rsub|h>-P<rsub|H>><rsub|1><l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|0>>>|<row|<cell|\<Rightarrow\><l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|0>>|<cell|\<leqslant\>>|<cell|<wide*|C<rsub|3>|\<wide-underbrace\>><rsub|=C<rsub|1>C<rsub|2>C<rsub|Q>>h<l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>>.>>>>
    </eqnarray*>

    <em|2nd step.>

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<lambda\><rsub|h>>|<cell|=>|<cell|sup<rsub|<l2norm|v<rsub|h>><rsub|0>=1><l2norm|A<rsub|h>v<rsub|h>><rsub|0>=sup<rsub|<l2norm|v<rsub|h>><rsub|0>=1>sup<rsub|<l2norm|w<rsub|h>><rsub|0>=1><wide*|<ip|A<rsub|h>v<rsub|h>|w<rsub|h>>|\<wide-underbrace\>><rsub|=a(v<rsub|h>,w<rsub|h>)>\<leqslant\>C<rsub|1>sup<rsub|<l2norm|v<rsub|h>><rsub|0>=<l2norm|w<rsub|h>><rsub|0>=1><l2norm|v<rsub|h>><rsub|1><l2norm|w<rsub|h>><rsub|1>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|C<rsub|1>C<rsub|I><rsup|2>h<rsup|-2>\<leqslant\>C<rsub|I>h<rsup|-1><l2norm|w<rsub|h>><rsub|0>>>|<row|<cell|\<Rightarrow\>h>|<cell|\<leqslant\>>|<cell|<frac|1|<sqrt|C<rsub|1>C<rsub|I><rsup|2>>>\<lambda\><rsub|h><rsup|-1/2>>>>>
    </eqnarray*>

    <with|color|red|Fehlt hier was?>

    <em|3rd step.>

    <\eqnarray*>
      <tformat|<table|<row|<cell|c<rsub|0><l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|1><rsup|2>>|<cell|\<leqslant\>>|<cell|a(v<rsub|h>-P<rsub|H>v<rsub|h>,v<rsub|h>-P<rsub|H>v<rsub|h>)>>|<row|<cell|>|<cell|=>|<cell|a(v<rsub|h>-P<rsub|H>v<rsub|h>,v<rsub|h>)>>|<row|<cell|>|<cell|=>|<cell|<ip|A<rsub|h>(v<rsub|h>-P<rsub|H>v<rsub|h>)|v<rsub|h>>>>|<row|<cell|>|<cell|=>|<cell|<ip|v<rsub|h>-P<rsub|H>|A<rsup|T><rsub|h>v<rsub|h>><rsub|0>\<leqslant\><l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|0><wide*|<l2norm|A<rsup|T><rsub|h>v<rsub|h>>|\<wide-underbrace\>><rsub|\<leqslant\>C<l2norm|A<rsub|h>v<rsub|h>>>\<leqslant\>C<rsub|4><rsub|>\<lambda\><rsup|-1/2><rsub|h><l2norm|A<rsub|h>v<rsub|h>><l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|1>>>>>
    </eqnarray*>

    \;
  </proof>

  <\remark*>
    Consider example <inactive|<reference|ex:0>> with <with|mode|math|h=1/N>,
    <with|mode|math|X<rsub|h>={v\<in\>C[0,1]:v<with|mode|text| linear in
    <with|mode|math|[x<rsub|n-1>,x<rsub|n>]>>}.> For
    <with|mode|math|v\<in\>C<rsup|2>[0,1]> define
    <with|mode|math|v<rsub|h>\<in\>X<rsub|h>> by
    <with|mode|math|v<rsub|h>(x<rsub|h>)=v(x<rsub|h>)>,
    <with|mode|math|x<rsub|n>=n*h>, show <with|mode|math|<wide*|<l2norm|v-v<rsub|h>><rsub|1>|\<wide-underbrace\>><rsub|=<l2norm|v<rprime|'>-v<rprime|'><rsub|h>><rsub|0>>\<leqslant\>C*h<l2norm|v<rprime|''>><rsub|0>>
    . Set <with|mode|math|w\<assign\>v-v<rsub|h>>,
    <with|mode|math|w(x<rsub|n>)=0> (<with|mode|math|n=0,\<ldots\>,N)>. There
    exists <with|mode|math|\<xi\><rsub|n>\<in\>[x<rsub|n-1>,x<rsub|n>]> with
    <with|mode|math|w<rprime|'>(\<xi\><rsub|n>)=0>.

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<Rightarrow\><l2norm|v<rprime|'>-v<rsub|h><rprime|'>><rsub|0><rsup|2>>|<cell|=>|<cell|<l2norm|w<rprime|'>><rsub|0><rsup|2>=<big|int><rsub|0><rsup|1>(w<rprime|'>)<rsup|2>d
      t=<big|sum><rsub|n=1><rsup|N><big|int><rsub|x<rsub|n-1>><rsup|x<rsub|n>><left|(>w<rprime|'>(t)-w<rprime|'>(\<xi\><rsub|n>)<right|)><rsup|2>d
      t>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<big|sum><rsub|n=1><rsup|N><big|int><rsub|x<rsub|n-1>><rsup|x<rsub|n>><left|(><big|int><rsub|\<xi\><rsub|n>><rsup|t>w<rprime|''>(s)d
      s<right|)><rsup|2>d t>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<big|sum><rsub|n=1><rsup|N><big|int><rsub|x<rsub|n-1>><rsup|x<rsub|n>><big|int><rsub|\<xi\><rsub|n>><rsup|t>d
      s<big|int><rsub|\<xi\><rsub|n>><rsup|t>(w<rprime|''>(s))<rsup|2>d s*d
      t>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<big|sum><rsub|n=1><rsup|N><big|int><rsub|x<rsub|n-1>><rsup|x<rsub|n>>\|t-\<xi\><rsub|n>\|<big|int><rsub|x<rsub|n-1>><rsup|x<rsub|n>>(w<rprime|''>(s))<rsup|2>d
      s*d t>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<frac|h|2><l2norm|w<rprime|''>><rsub|0><rsup|2>>>>>
    </eqnarray*>

    ad b)

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|v<rsub|h><rprime|'>><rsub|0><rsup|2>>|<cell|=>|<cell|<big|sum><rsub|n=1><rsup|N><big|int><rsub|x<rsub|n-1>><rsup|x<rsub|n>><left|(><frac|1|h>(v(x<rsub|n>)-v(x<rsub|n-1>))<right|)><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|h<rsup|-1><big|sum><rsub|n=1><rsup|N>(v(x<rsub|n>)-v(x<rsub|n-1>))<rsup|2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|2h<rsup|-2><big|sum><rsub|n=1><rsup|N>(v(x<rsub|n>)+v(x<rsub|n-1>))<rsup|2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<frac|12|h<rsup|2>><big|sum><rsub|n=1><rsup|N>h<frac|1|6><left|(>v(x<rsub|n-1>)<rsup|2>+4<left|(><frac|v(x<rsub|n>)+v(x<rsub|n-1>)|2><right|)><rsup|2>+v(x<rsub|n>)<rsup|2><right|)>=<l2norm|v<rsub|>><rsup|2><rsub|0>.>>>>
    </eqnarray*>
  </remark*>

  <\corollary>
    For symmetric <with|mode|math|A> and <with|mode|math|<l2norm|\<cdot\>><rsub|1>=<enorm|\<cdot\>>>
    the following are equivalent:

    <\enumerate-numeric>
      <item><with|mode|math|<l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|0>\<leqslant\>C*H<l2norm|v<rsub|h>-P<rsub|H>><rsub|1>>,

      <item><with|mode|math|<l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|1>\<leqslant\>C*H<l2norm|A<rsub|h>v<rsub|h>><rsub|0>>,

      <item>if <with|mode|math|<l2norm|v<rsub|h>><rsub|1>\<leqslant\>C<rsub|I>h<rsup|-1><l2norm|v<rsub|h>><rsub|0>>
      and <with|mode|math|<l2norm|v<rsub|h>-Q<rsub|H>v<rsub|h>><rsub|0>=C*H<l2norm|v<rsub|h>>>,
      <with|mode|math|<l2norm|P<rsub|H>v<rsub|h>><rsub|0>\<leqslant\>C<l2norm|v<rsub|h>><rsub|0>>.
    </enumerate-numeric>
  </corollary>

  <\proof>
    1<with|mode|math|\<Rightarrow\>>2: done.

    2<with|mode|math|\<Rightarrow\>>1:

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>>>|<cell|=>|<cell|sup<rsub|w<rsub|h>\<neq\>0><frac|<ip|v<rsub|h>-P<rsub|H>v<rsub|h>|w<rsub|h>>|<l2norm|w<rsub|h>><rsub|0>>=sup<rsub|w<rsub|h>\<neq\>0><frac|<ip|v<rsub|h>-P<rsub|H>v<rsub|h>|A*w<rsub|h>>|<l2norm|A*w<rsub|h>><rsub|0>>>>|<row|<cell|>|<cell|=>|<cell|sup<rsub|w<rsub|h>\<neq\>0><frac|a(v<rsub|h>-P<rsub|H>v<rsub|h>,w<rsub|h>)|<l2norm|A*w<rsub|h>><rsub|0>>\<leqslant\>sup<rsub|w<rsub|h>\<neq\>0><frac|a(v<rsub|h>-P<rsub|H>v<rsub|h>,w<rsub|h>-P<rsub|H>w<rsub|h>)|<l2norm|A*w<rsub|h>><rsub|0>>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|sup<rsub|w<rsub|h>\<neq\>0>C<rsub|1><l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|1><wide*|<frac|<l2norm|w<rsub|h>-P<rsub|H>w<rsub|h>>|<l2norm|A<rsub|h>w<rsub|h>>>|\<wide-underbrace\>><rsub|=C*H>.>>>>
    </eqnarray*>

    1<with|mode|math|\<Rightarrow\>3:>

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|P<rsub|H>v<rsub|h>><rsub|0>>|<cell|\<leqslant\>>|<cell|<l2norm|v<rsub|h>><rsub|0>+<l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|0>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<l2norm|v<rsub|h>><rsub|0>+C*H<l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|1>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<l2norm|v<rsub|h>><rsub|0>+C*H<l2norm|v<rsub|h>><rsub|1>\<leqslant\>C*HC<rsub|I>h<rsup|-1><l2norm|v<rsub|h>><rsub|0>.>>>>
    </eqnarray*>

    3<with|mode|math|\<Rightarrow\>>1:

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|A<rsub|h>P<rsub|H>v<rsub|h>><rsub|0><rsup|2>>|<cell|=>|<cell|<ip|A<rsub|h>P<rsub|H>v<rsub|h>|A<rsub|h>P<rsub|H>v<rsub|h>>=a(P<rsub|H>v<rsub|h>,A<rsub|h>P<rsub|H>v<rsub|h>)>>|<row|<cell|>|<cell|=>|<cell|a(v<rsub|h>,P<rsub|H>A<rsub|h>P<rsub|H>v<rsub|h>)=<ip|A<rsub|h>v<rsub|h>|P<rsub|H>A<rsub|h>P<rsub|H>v<rsub|h>>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<l2norm|A<rsub|h>v<rsub|h>><rsub|0><l2norm|P<rsub|H>A<rsub|h>P<rsub|H>v<rsub|h>><rsub|0>\<leqslant\>C<l2norm|A<rsub|h>P<rsub|H>v<rsub|h>><rsub|0>.>>>>
    </eqnarray*>

    Now,

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|0>>|<cell|=>|<cell|sup<rsub|w<rsub|h>\<neq\>0><frac|<ip|v<rsub|h>-P<rsub|H>v<rsub|h>|A*w<rsub|h>>|<l2norm|A<rsub|h>w<rsub|h>>>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|sup<rsub|w<rsub|h>\<neq\>0><frac|C*H<l2norm|v<rsub|h>-P<rsub|H>v<rsub|h>><rsub|0>(1+C)<l2norm|A<rsub|h>w<rsub|h>><rsub|0>|<l2norm|A<rsub|h>w<rsub|h>>>,>>>>
    </eqnarray*>

    considering

    <\eqnarray*>
      <tformat|<table|<row|<cell|a(v<rsub|h>-P<rsub|H>v<rsub|h>,w<rsub|h>)>|<cell|=>|<cell|a(v<rsub|h>-P<rsub|H>v<rsub|h>,w<rsub|h>-P<rsub|H>w<rsub|h>)>>|<row|<cell|>|<cell|=>|<cell|a((v<rsub|h>-P<rsub|H>v<rsub|h>)-Q<rsub|H>(v<rsub|h>-P<rsub|H>v<rsub|h>),w<rsub|h>-P<rsub|H>w<rsub|h>)>>|<row|<cell|>|<cell|=>|<cell|<ip|(v<rsub|h>-P<rsub|H>v<rsub|h>)-Q<rsub|H>(v<rsub|h>-P<rsub|H>v<rsub|h>)|A<rsub|h>(w<rsub|h>-P<rsub|H>w<rsub|h>)><rsub|0>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<l2norm|(v<rsub|h>-P<rsub|H>v<rsub|h>)-Q<rsub|H>(v<rsub|h>-P<rsub|H>v<rsub|h>)><rsub|0><left|(><l2norm|A<rsub|h>w<rsub|h>><rsub|0>+<l2norm|A<rsub|h>P<rsub|H>w<rsub|h>><rsub|0><right|)>>>>>
    </eqnarray*>

    \;
  </proof>

  <\lemma>
    Let <with|mode|math|A<rsub|h>> and <with|mode|math|\<theta\><rsub|h>\<leqslant\>\<lambda\><rsub|h><rsup|-1>>.
    Then, the smoothing property

    <\equation*>
      <l2norm|A<rsub|h>(id-\<theta\><rsub|h>A<rsub|h>)<rsup|m>>\<leqslant\>\<lambda\><rsub|h><frac|C<rsub|R>|m>
    </equation*>

    holds.
  </lemma>

  <\proof>
    <with|mode|math|A<rsub|h>w<rsub|h><rsup|j>=\<mu\><rsub|j>w<rsub|h><rsup|j>>,
    with a complete ONS <with|mode|math|w<rsub|h><rsup|j>>.

    <\equation*>
      v<rsub|h>=<big|sum><ip|w<rsub|h><rsup|j>|v<rsub|h>>w<rsub|h><rsup|j>.
    </equation*>

    Then

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|A<rsub|h>(id-\<theta\><rsub|h>A<rsub|h>)<rsup|m>v<rsub|h>><rsup|2>>|<cell|=>|<cell|<l2norm|<big|sum><ip|w<rsub|h><rsup|j>|v<rsub|h>>A<rsub|h>(id-\<theta\><rsub|h>A<rsub|h>)<rsup|m>w<rsub|h><rsup|j>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|\<theta\><rsub|h><rsup|-1><l2norm|<big|sum><ip|w<rsub|h><rsup|j>|v<rsub|h>><rsub|0><wide*|\<theta\><rsub|h>\<mu\><rsub|j>|\<wide-underbrace\>><rsub|t>(1-\<theta\><rsub|h>\<mu\><rsub|j>)<rsup|m>w<rsub|j>><rsub|0>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|\<ldots\>\<ldots\>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|\<lambda\><rsub|h>max<rsub|t\<in\>[0,1]><left|(>t(1-t)<rsup|m><right|)><rsup|2><l2norm|v<rsub|h>><rsub|0><rsup|2>,>>>>
    </eqnarray*>

    considering <with|mode|math|\<lambda\><rsub|h>=max \<mu\><rsub|j>>.
  </proof>
</body>

<\initial>
  <\collection>
    <associate|language|american>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|3|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Classical
      Two-Level Analysis> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>