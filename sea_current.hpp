#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <numbers>
#include <algorithm>
#include <tuple>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>


#ifdef DEBUG
inline void dbg_assert(bool cnd, const char* msg) {
    if (!cnd) {
        std::cerr << "assertion failed: " << msg << std::endl;
        std::exit(1);
    }
}
#else
#define dbg_assert(a, b)
#endif

namespace turtle::sc {

    using namespace Eigen;

    using namespace std::complex_literals;


    struct arclength_data {
        float arclength;
        std::vector<VectorXf> segments;
        std::vector<VectorXf> positions;
    };

    struct chebpoly {
        VectorXf coeffs;
        float xmin;
        float xmax;
        chebpoly(VectorXf coeffs, const float xmin, const float xmax) : coeffs(coeffs), xmin(xmin), xmax(xmax) {}
    };

    chebpoly chebfit(const VectorXf& x, const VectorXf& y, const int degree);
    VectorXf chebeval(const VectorXf& x, const chebpoly& b, const int degree);

    class bezier_spline {
        public:
            std::vector<std::vector<Vector2f>> ctrl_pts;
            Matrix<float, Dynamic, 2> pts;
            std::vector<VectorXf> positions; // for curvature/hodograph

            bezier_spline() = default;
            bezier_spline(const std::vector<std::vector<Vector2f>>& ctrl_pts, const Matrix<float, Dynamic, 2>& pts, const std::vector<VectorXf>& positions);

            static bezier_spline join_splines(std::vector<bezier_spline>& splines);

            static bezier_spline bezier_curve(const std::vector<Vector2f>& ctrl_pts, const std::vector<float>& positions);
            static bezier_spline bezier_curve(const std::vector<Vector2f>& ctrl_pts, const VectorXf& positions);
            static bezier_spline bezier_curve(std::vector<Vector2f>& ctrl_pts, float precision);

            inline int n_pts() const;
            inline int n_segments() const;
            inline int degree() const;

            arclength_data arclength(const float precision) const;
            std::vector<float> curvature() const;
            bezier_spline hodograph() const;
            bezier_spline resample(VectorXf& profile_pos, arclength_data ad, bool nudge_positions) const;

            // helper function for bezier_curve
            static std::vector<std::complex<float>> omega_table(int degree);

        private:
            template <typename T, typename Y>
            static inline T coerce(T& num, Y low, Y high) {
                if (num < low) return low;
                else if (num > high) return high;
                return num;
            }
    };

    bezier_spline::bezier_spline(const std::vector<std::vector<Vector2f>>& ctrl_pts, const Matrix<float, Dynamic, 2>& pts, const std::vector<VectorXf>& positions) : ctrl_pts(ctrl_pts), pts(pts), positions(positions) {}

    inline int bezier_spline::n_pts() const {
        return pts.rows();
    }

    inline int bezier_spline::n_segments() const {
        return ctrl_pts.size();
    }

    inline int bezier_spline::degree() const {
        return ctrl_pts[0].size()-1;
    }

    // it's ok that this is static as a new matrix would have to be allocated either way
    bezier_spline bezier_spline::join_splines(std::vector<bezier_spline>& splines) {
        std::vector<std::vector<Vector2f>> ctrl_pts;
        std::vector<VectorXf> positions;

        int n_pts = 0;
        for (bezier_spline spline : splines) {
            n_pts += spline.n_pts();
            for (std::size_t i = 0; i < spline.n_segments(); ++i) {
                ctrl_pts.push_back(spline.ctrl_pts[i]);
                positions.push_back(spline.positions[i]);
            }
        }

        Matrix<float, Dynamic, 2> joined_pts = Matrix<float, Dynamic, 2>::Zero(n_pts, 2);

        int n = 0;
        for (bezier_spline spline : splines) {
            joined_pts.block(n, 0, spline.n_pts(), 2) = spline.pts;
            n += spline.n_pts();
        }

        return bezier_spline(ctrl_pts, joined_pts, positions);
    }

    bezier_spline bezier_spline::bezier_curve(const std::vector<Vector2f>& ctrl_pts, const std::vector<float>& positions) {
        const VectorXf mapped = VectorXf::Map(&positions[0], positions.size());
        return bezier_curve(ctrl_pts, mapped);
    }

    bezier_spline bezier_spline::bezier_curve(const std::vector<Vector2f>& ctrl_pts, const VectorXf& positions) {
        dbg_assert(ctrl_pts.size() >= 2, "ctrl_pts must have at least 2 points");
        dbg_assert(positions.minCoeff() >= 0, "positions must be between [0, 1]");
        dbg_assert(positions.minCoeff() <= 1, "positions must be between [0, 1]");

        const int degree = ctrl_pts.size() - 1;

        const std::vector<std::complex<float>> omegas = bezier_spline::omega_table(degree);

        FFT<float> fft;

        Matrix<std::complex<float>, Dynamic, 2> U = Matrix<float, Dynamic, 2>::Zero(ctrl_pts.size(), 2);

        for (std::size_t i = 0; i <= degree; ++i) {
            U(i, 0) = ctrl_pts[i].x();
            U(i, 1) = ctrl_pts[i].y();
        }

        Matrix<std::complex<float>, Dynamic, 2> Q = Matrix<float, Dynamic, 2>::Zero(ctrl_pts.size(), 2);

        Eigen::VectorXcf tmp_fft(ctrl_pts.size());
        fft.inv(tmp_fft, U.col(0));
        Q.col(0) = tmp_fft;

        fft.inv(tmp_fft, U.col(1));
        Q.col(1) = tmp_fft;

        Matrix<float, Dynamic, 2> B;
        B.setZero(positions.rows(), 2);

        for (int i = 0; i < positions.rows(); ++i) {
            const float s = positions(i);
            for (int k = 0; k <= degree; ++k) {
                std::complex<float> tmp = std::pow((1.0f+0if) + s*(omegas[k] - (1.0f+0if)), degree);
                B(i, 0) += (Q(k, 0) * tmp).real();
                B(i, 1) += (Q(k, 1) * tmp).real();
            }
        }

        return bezier_spline({ctrl_pts}, B, {positions});
    }

    bezier_spline bezier_spline::bezier_curve(std::vector<Vector2f>& ctrl_pts, const float precision) {
        dbg_assert(precision < 1 && precision > 0, "spline percision must be in (0, 1)");

        const int n = static_cast<int>(std::round(1.0f / precision));

        dbg_assert(std::abs(n - (1.0f / precision)) < 0.0001,
                    "1/percision must be (close) to an integer, for arbitrary position values use the other bezier_curve method");

        std::vector<float> positions(n+1);
        for (std::size_t i = 0; i <= n; ++i) {
            positions[i] = std::min(i * precision, 1.0f);
        }

        return bezier_curve(ctrl_pts, positions);
    }


    // std::tuple<float, std::vector<std::vector<float>>> bezier_spline::arclength(const float precision=0.01) const {
    arclength_data bezier_spline::arclength(const float precision=0.01) const {
        constexpr std::array<double, 32> weights = {
            0.0965400885147278005667648300635757947368606312355700687323182099577497758679466512968173871061464644599963197828969869820251559172455698832434930732077927850876632725829187045819145660710266452161095406358159608874152584850413283587913891015545638518881205600825069096855488296437485836866,
            0.0965400885147278005667648300635757947368606312355700687323182099577497758679466512968173871061464644599963197828969869820251559172455698832434930732077927850876632725829187045819145660710266452161095406358159608874152584850413283587913891015545638518881205600825069096855488296437485836866,
            0.0956387200792748594190820022041311005948905081620055509529898509437067444366006256133614167190847508238474888230077112990752876436158047205555474265705582078453283640212465537132165041268773645168746774530146140911679782502276289938840330631903789120176765314495900053061764438990021439069,
            0.0956387200792748594190820022041311005948905081620055509529898509437067444366006256133614167190847508238474888230077112990752876436158047205555474265705582078453283640212465537132165041268773645168746774530146140911679782502276289938840330631903789120176765314495900053061764438990021439069,
            0.0938443990808045656391802376681172600361000757462364500506275696355695118623098075097804207682530277555307864917078828352419853248607668520631751470962234105835015158485760721979732297206950719908744248285672032436598213262204039212897239890934116841559005147755270269705682414708355646603,
            0.0938443990808045656391802376681172600361000757462364500506275696355695118623098075097804207682530277555307864917078828352419853248607668520631751470962234105835015158485760721979732297206950719908744248285672032436598213262204039212897239890934116841559005147755270269705682414708355646603,
            0.0911738786957638847128685771116370625448614132753900053231278739777031520613017513597426417145878622654027367650308019870251963114683369110451524174258161390823876554910693202594383388549640738095422966058367070348943662290656339592299608558384147559830707904449930677260444604329157917977,
            0.0911738786957638847128685771116370625448614132753900053231278739777031520613017513597426417145878622654027367650308019870251963114683369110451524174258161390823876554910693202594383388549640738095422966058367070348943662290656339592299608558384147559830707904449930677260444604329157917977,
            0.0876520930044038111427714627518022875484497217017572223192228034747061150211380239263021665771581379364685191248848158059408000065275041643745927401342920150588893827207354226012701872322225514682178439577327346929209121046816487338309068375228210705166692551938339727096609740531893725675,
            0.0876520930044038111427714627518022875484497217017572223192228034747061150211380239263021665771581379364685191248848158059408000065275041643745927401342920150588893827207354226012701872322225514682178439577327346929209121046816487338309068375228210705166692551938339727096609740531893725675,
            0.0833119242269467552221990746043486115387468839428344598401864047287594069244380966536255650452315042012372905572506028852130723585016898197140339352228963465326746426938359210160503509807644396182380868089959855742801355208471205261406307895519604387550841954817025499019984032594036141439,
            0.0833119242269467552221990746043486115387468839428344598401864047287594069244380966536255650452315042012372905572506028852130723585016898197140339352228963465326746426938359210160503509807644396182380868089959855742801355208471205261406307895519604387550841954817025499019984032594036141439,
            0.078193895787070306471740918828306671039786798482159190307481553869493700115196435401943819761440851294456424770323467367505109006517482028994114252939401250416132320553639542341400437522236191275346323130525969269563653003188829786549728825182082678498917784036375053244425839341945385297,
            0.078193895787070306471740918828306671039786798482159190307481553869493700115196435401943819761440851294456424770323467367505109006517482028994114252939401250416132320553639542341400437522236191275346323130525969269563653003188829786549728825182082678498917784036375053244425839341945385297,
            0.0723457941088485062253993564784877916043369833018248707397632823511765345816800402874475958591657429073027694582930574378890633404841054620298756279975430795706338162404545590689277985270140590721779502609564199074051863640176937117952488466002340085264819537808079947788437998042296495822,
            0.0723457941088485062253993564784877916043369833018248707397632823511765345816800402874475958591657429073027694582930574378890633404841054620298756279975430795706338162404545590689277985270140590721779502609564199074051863640176937117952488466002340085264819537808079947788437998042296495822,
            0.0658222227763618468376500637069387728775364473732465153710916696852412442018627316280044447764609054151761388378861151807154113495715653711918644796313239555117970398473141615070299152284100887258072240524028885129828725430021172354299810423059697133688823072212214503334259555369485963074,
            0.0658222227763618468376500637069387728775364473732465153710916696852412442018627316280044447764609054151761388378861151807154113495715653711918644796313239555117970398473141615070299152284100887258072240524028885129828725430021172354299810423059697133688823072212214503334259555369485963074,
            0.0586840934785355471452836373001708867501204674575467587150032786132877518019090643743123653437052116901895704813134467814193905269714480573030647540887991405215103758723074481312705449946311993670933802369300463315125015975216910705047901943865293781921122370996257470349807212516159332678,
            0.0586840934785355471452836373001708867501204674575467587150032786132877518019090643743123653437052116901895704813134467814193905269714480573030647540887991405215103758723074481312705449946311993670933802369300463315125015975216910705047901943865293781921122370996257470349807212516159332678,
            0.0509980592623761761961632446895216952601847767397628437069071236525030510385137821267442193868358292147899714519363571211100873456269865150186456681043804358654826791768545393024953758025593924464295555854744882720755747096079325496814455853004350452095212995888025282619932613606999567133,
            0.0509980592623761761961632446895216952601847767397628437069071236525030510385137821267442193868358292147899714519363571211100873456269865150186456681043804358654826791768545393024953758025593924464295555854744882720755747096079325496814455853004350452095212995888025282619932613606999567133,
            0.0428358980222266806568786466061255284928108575989407395620219408911043916962572261359138025961596979511472539467367407419206021900868371610612953162236233351132214438513203223655531564777278515080476421262443325932320214191168239648611793958596884827086182431203349730049744697408543115307,
            0.0428358980222266806568786466061255284928108575989407395620219408911043916962572261359138025961596979511472539467367407419206021900868371610612953162236233351132214438513203223655531564777278515080476421262443325932320214191168239648611793958596884827086182431203349730049744697408543115307,
            0.0342738629130214331026877322523727069948402029116274337814057454192310522168984446294442724624445760666244242305266023810860790282088335398182296698622433517061843276344829146573593201201081743714879684153735672789104567624853712011151505225193933019375481618760594889854480408562043658635,
            0.0342738629130214331026877322523727069948402029116274337814057454192310522168984446294442724624445760666244242305266023810860790282088335398182296698622433517061843276344829146573593201201081743714879684153735672789104567624853712011151505225193933019375481618760594889854480408562043658635,
            0.0253920653092620594557525897892240292875540475469487209362512822192154788532376645960457016338988332029324531233401833547954942765653767672102838323550828207273795044402516181251040411735351747299230615776597356956641506445501689924551185923348003766988424170511157069264716719906995309826,
            0.0253920653092620594557525897892240292875540475469487209362512822192154788532376645960457016338988332029324531233401833547954942765653767672102838323550828207273795044402516181251040411735351747299230615776597356956641506445501689924551185923348003766988424170511157069264716719906995309826,
            0.0162743947309056706051705622063866181795429637952095664295931749613369651752917857651844425586692833071042366002861684552859449530958901379260437604156888337987656773068694383447504913457771896770689760342192010638946676879735404121702279005140285599424477022083127753774756520463311689155,
            0.0162743947309056706051705622063866181795429637952095664295931749613369651752917857651844425586692833071042366002861684552859449530958901379260437604156888337987656773068694383447504913457771896770689760342192010638946676879735404121702279005140285599424477022083127753774756520463311689155,
            0.0070186100094700966004070637388531825133772207289396032320082356192151241454178686953297376907573215077936155545790593837513204206518026084505878987243348925784479817181234617862457418214505322067610482902501455504204433524520665822704844582452877416001060465891907497519632353148380799619,
            0.0070186100094700966004070637388531825133772207289396032320082356192151241454178686953297376907573215077936155545790593837513204206518026084505878987243348925784479817181234617862457418214505322067610482902501455504204433524520665822704844582452877416001060465891907497519632353148380799619};

        constexpr std::array<double, 32> abscissa = {
            -0.0483076656877383162348125704405021636908472517308488971677937345463685926042778777794060365911173780988289503411375793689757446357461295741679964108035347980667582792392651327368009453047606446744575790523465655622949909588624860214137051585425884056992683442137333250625173849291299678673,
            0.0483076656877383162348125704405021636908472517308488971677937345463685926042778777794060365911173780988289503411375793689757446357461295741679964108035347980667582792392651327368009453047606446744575790523465655622949909588624860214137051585425884056992683442137333250625173849291299678673,
            -0.1444719615827964934851863735988106522038459913156355521379528938242184438164519731102406769974924713989580220758441301598578946580142268413547299935841673092513202403499286272686350814272974392746706128556678811982653393383080797337231702069432462445053984587997153683967433095128570624414,
            0.1444719615827964934851863735988106522038459913156355521379528938242184438164519731102406769974924713989580220758441301598578946580142268413547299935841673092513202403499286272686350814272974392746706128556678811982653393383080797337231702069432462445053984587997153683967433095128570624414,
            -0.2392873622521370745446032091655015206088554219602530155470960995597029133039943915553593695844147813728958071901224632260145752503694970545640339873418480550362677768010887468668377893757173424222709744116861683634989914911762187599464033126988486345234374380695224452457957624756811128321,
            0.2392873622521370745446032091655015206088554219602530155470960995597029133039943915553593695844147813728958071901224632260145752503694970545640339873418480550362677768010887468668377893757173424222709744116861683634989914911762187599464033126988486345234374380695224452457957624756811128321,
            -0.3318686022821276497799168057301879961957751368050598360182296306285376829657438169809731852312743263005943551508559377834274303920771100489026913715847854727626540340157368609696698131829681988642689780208633461925468064919389286805624602715005948661328152252049795463242055567997437182143,
            0.3318686022821276497799168057301879961957751368050598360182296306285376829657438169809731852312743263005943551508559377834274303920771100489026913715847854727626540340157368609696698131829681988642689780208633461925468064919389286805624602715005948661328152252049795463242055567997437182143,
            -0.4213512761306353453641194361724264783358772886324433305416613404557190462549837315607633055675740638739884093394574651160978879545562247406839036854173715776910866941643197988581928900702286425821151586000969947406313405310082646561917980302543820974679501841964453794193724645925031841919,
            0.4213512761306353453641194361724264783358772886324433305416613404557190462549837315607633055675740638739884093394574651160978879545562247406839036854173715776910866941643197988581928900702286425821151586000969947406313405310082646561917980302543820974679501841964453794193724645925031841919,
            -0.5068999089322293900237474743778212301802836995994354639743662809707712640478764442266190213124522047999876916596854537447047905434649918210338296049592120273725464263651562560829050004258268002241145951271730860506703690843719936432852920782304931272053564539127514959875734718036950073563,
            0.5068999089322293900237474743778212301802836995994354639743662809707712640478764442266190213124522047999876916596854537447047905434649918210338296049592120273725464263651562560829050004258268002241145951271730860506703690843719936432852920782304931272053564539127514959875734718036950073563,
            -0.5877157572407623290407454764018268584509401154544205727031788473129228586684474311408145102018661764979429510790747919023774933113319119601088669936958908618326367715806216053155906936017362413244183150445492317940727345571648726363597097311647731726438279098059670236086983675374932643925,
            0.5877157572407623290407454764018268584509401154544205727031788473129228586684474311408145102018661764979429510790747919023774933113319119601088669936958908618326367715806216053155906936017362413244183150445492317940727345571648726363597097311647731726438279098059670236086983675374932643925,
            -0.6630442669302152009751151686632383689770222859605053010170834964924461749232229404368981536611965356686820332804126742949900731319113817214392193185613161549689934301410316417342588149871686184296988807305719690974644891055567340650986465615021143958920599684258616066247948224049997371166,
            0.6630442669302152009751151686632383689770222859605053010170834964924461749232229404368981536611965356686820332804126742949900731319113817214392193185613161549689934301410316417342588149871686184296988807305719690974644891055567340650986465615021143958920599684258616066247948224049997371166,
            -0.732182118740289680387426665091267146630270483506629100821139573270385253587797727611292298988652560055905228466313310601075333829094630570926240639601009902567982815376254840388565733846030450161774620971196087756484387383432502715118096615117242484073636640563609696801484680439912327302,
            0.732182118740289680387426665091267146630270483506629100821139573270385253587797727611292298988652560055905228466313310601075333829094630570926240639601009902567982815376254840388565733846030450161774620971196087756484387383432502715118096615117242484073636640563609696801484680439912327302,
            -0.7944837959679424069630972989704289020954794016388354532507582449720593922816426654241878967890821228397041480126630294067578180914548706957761322921470535094589673860419616615738928385807346185892317514562489971543238450942224396667500582904031225063621511429185567036727089257387570529468,
            0.7944837959679424069630972989704289020954794016388354532507582449720593922816426654241878967890821228397041480126630294067578180914548706957761322921470535094589673860419616615738928385807346185892317514562489971543238450942224396667500582904031225063621511429185567036727089257387570529468,
            -0.849367613732569970133693004967742538954886793049759233100219598613724656141562558741881463752754991143937635778596582088915769685796612254240615386941355933272723068952531445772190363422003834495043219316062885999846179078139659341918527603834809670576387535564876596379488780285979062125,
            0.849367613732569970133693004967742538954886793049759233100219598613724656141562558741881463752754991143937635778596582088915769685796612254240615386941355933272723068952531445772190363422003834495043219316062885999846179078139659341918527603834809670576387535564876596379488780285979062125,
            -0.8963211557660521239653072437192122684789964967957595765636154129650249794910409173494503783167666654202705333374285522819507600044591355080910768854012859468015827508424619812224062460791781333400979810176198916239783226706506012473250929962326307746466256167673927887144428859779028909399,
            0.8963211557660521239653072437192122684789964967957595765636154129650249794910409173494503783167666654202705333374285522819507600044591355080910768854012859468015827508424619812224062460791781333400979810176198916239783226706506012473250929962326307746466256167673927887144428859779028909399,
            -0.9349060759377396891709191348354093255286714322828372184584037398118161947182932855418880831417927728359606280450921427988850058691931014887248988124656348299653052688344696135840215712191162135178273756415771123010111796122671724143565383396162107206772781551029308751511942924942333859805,
            0.9349060759377396891709191348354093255286714322828372184584037398118161947182932855418880831417927728359606280450921427988850058691931014887248988124656348299653052688344696135840215712191162135178273756415771123010111796122671724143565383396162107206772781551029308751511942924942333859805,
            -0.9647622555875064307738119281182749603888952204430187193220113218370995254867038008243801877562227002840740910741483519987441236283464394249183812395373150090695515823078220949436846111682404866338388944248976976566275875721000356873959697266702651250019105084704924793016185368873243713355,
            0.9647622555875064307738119281182749603888952204430187193220113218370995254867038008243801877562227002840740910741483519987441236283464394249183812395373150090695515823078220949436846111682404866338388944248976976566275875721000356873959697266702651250019105084704924793016185368873243713355,
            -0.9856115115452683354001750446309019786323957143358063182107821705820305847193755946663846485510970266115353839862364606643634021712823093784875255943834038377710426488328772047833289470320023596895438028281274741367781028592272459887917924171204666683239464005128153533797603112851826904814,
            0.9856115115452683354001750446309019786323957143358063182107821705820305847193755946663846485510970266115353839862364606643634021712823093784875255943834038377710426488328772047833289470320023596895438028281274741367781028592272459887917924171204666683239464005128153533797603112851826904814,
            -0.9972638618494815635449811286650407271385376637294611593011185457862359083917418520130456693085426416474280482200936551645510686196373231416035137741332968299789863385253514914078766236061488136738023162574655835389902337937054326098485227311719825229066712510246574949376367552421728646398,
            0.9972638618494815635449811286650407271385376637294611593011185457862359083917418520130456693085426416474280482200936551645510686196373231416035137741332968299789863385253514914078766236061488136738023162574655835389902337937054326098485227311719825229066712510246574949376367552421728646398};

        static_assert(weights.size() == abscissa.size());

        float total_arclen = 0;
        std::vector<VectorXf> arclens(n_segments());
        std::vector<VectorXf> arclen_positions(n_segments());
        for (int s = 0; s < n_segments(); ++s) {
            dbg_assert(precision < 1 && precision > 0, "spline percision must be in (0, 1)");

            const int n = static_cast<int>(std::round(1.0f / precision));

            dbg_assert(std::abs(n - (1.0f / precision)) < 0.0001,
                        "1/percision must be (close) to an integer, for arbitrary position values use the other bezier_curve method");
            
            std::vector<float> pos(n+1);
            for (std::size_t i = 0; i <= n; ++i) {
                pos[i] = std::min(i * precision, 1.0f);
            }

            VectorXf deriv_pos = VectorXf::Zero(weights.size() * (pos.size()-1));
            for (int k = 0; k < pos.size()-1; ++k) {
                const float b = pos[k+1];
                const float a = pos[k];

                std::vector<float> ab_pos(weights.size());
                for (int i = 0; i < weights.size(); ++i) {
                    deriv_pos((k*weights.size())+i) = ((b-a)/2)*abscissa[i] + ((b+a)/2);
                }
            }

            bezier_spline deriv = (bezier_spline({ctrl_pts[s]}, Matrix<float, Dynamic, 2>(), {deriv_pos})).hodograph();

            VectorXf arclen_seg = VectorXf::Zero(pos.size());
            VectorXf arclen_seg_pos = VectorXf::Zero(pos.size());
            for (int k = 0; k < pos.size()-1; ++k) {
                const float b = pos[k+1];
                const float a = pos[k];

                for (int i = 0; i < weights.size(); ++i) {
                    const int ind = (k*weights.size())+i;
                    arclen_seg(k+1) += weights[i] * std::sqrt((deriv.pts(ind, 0) * deriv.pts(ind, 0)) + (deriv.pts(ind, 1) * deriv.pts(ind, 1)));
                }
                arclen_seg[k+1] *= ((b-a)/2.0f);
                arclen_seg_pos[k+1] = b;
                total_arclen += arclen_seg[k+1];
            }

            float sum = 0;
            for (float& arc : arclen_seg) {
                sum += arc;
                arc = sum;
            }
            arclens[s] = arclen_seg;
            arclen_positions[s] = arclen_seg_pos;
        }
        arclength_data ad;
        ad.arclength = total_arclen;
        ad.segments = arclens;
        ad.positions = arclen_positions;

        return ad;
    }

    bezier_spline bezier_spline::resample(VectorXf& profile_pos, arclength_data ad, bool nudge_positions=false) const {

        // This method of nudging only works well for isolated cases of weird values
        if (nudge_positions) {
            profile_pos(0) = 0;
            profile_pos(profile_pos.rows()-1) = ad.arclength;
            for (int i = 1; i < profile_pos.rows()-1; ++i) {
                if (profile_pos(i) < profile_pos(i-1) || profile_pos(i) > profile_pos(i+1)) {
                    profile_pos(i) = (profile_pos(i-1) + profile_pos(i+1)) / 2;
                }
            }
        }

        dbg_assert(profile_pos.rows() > 0, "The vector of positions to be sampled must not be empty");
        dbg_assert(profile_pos.maxCoeff() <= ad.arclength, "The profile can not go beyond the arclength of the spline");
        dbg_assert(profile_pos.minCoeff() >= 0, "The profile can not go beyond the arclength of the spline");

        std::vector<bezier_spline> curves(positions.size());

        int j = 0;
        for (int i = 0; i < positions.size(); ++i) {
            const VectorXf seg = ad.segments[i];
            const float last = seg(seg.rows()-1);
            float offset = 0;
            int start = j;
            for (; j < profile_pos.rows() && (profile_pos(j) - offset <= last); ++j) {} // cursed
            j -= 1;
            offset = profile_pos(j);

            VectorXf block = profile_pos.block(start, 0, ((j+1)-start), 1).array() - profile_pos(start);

            const int degree = std::max(seg.rows()/2, 2L);


            // polynomial from segment arclength to [0,1]
            chebpoly poly = chebfit(seg, ad.positions[i], degree);
            VectorXf positions_fixed = chebeval(block, poly, degree);

            curves[i] = bezier_curve(ctrl_pts[i], positions_fixed);
        }

        bezier_spline fixed_spline = join_splines(curves);
        dbg_assert(fixed_spline.n_pts() == profile_pos.rows(), "fixed_spline.n_pts() == profile_pos.rows()");
        return fixed_spline;
    }

    std::vector<float> bezier_spline::curvature() const {


    }

    bezier_spline bezier_spline::hodograph() const {
        std::vector<bezier_spline> curves(ctrl_pts.size());

        for (std::size_t i = 0; i < ctrl_pts.size(); ++i) {
            std::vector<Vector2f> d_cps(ctrl_pts[i].size()-1);
            for (std::size_t j = 0; j < ctrl_pts[i].size()-1; ++j) {
                d_cps[j] = degree() * (ctrl_pts[i][j+1] - ctrl_pts[i][j]);
            }
            curves[i] = bezier_curve(d_cps, positions[i]);
        }

        return join_splines(curves);
    }

    std::vector<std::complex<float>> bezier_spline::omega_table(const int degree) {
        using std::complex;

        std::vector<std::complex<float>> omegas(degree+1);
        omegas[0] = 1.0f+0if;
        for (std::size_t i = 1; i <= degree; ++i) {
            omegas[i] = pow(exp((complex<float>(std::numbers::pi) * -2if) / complex<float>(degree+1)), i);
        }

        return omegas;
    }


    chebpoly chebfit(const VectorXf& x, const VectorXf& y, const int degree) {
        dbg_assert(degree >= 1, "degree must be a positive integer");
        dbg_assert(x.rows() == y.rows(), "x and y must have the same number of rows");

        const int n = degree;
        const int m = x.rows();
        const float xmax = x.maxCoeff();
        const float xmin = x.minCoeff();

        dbg_assert(std::abs(xmax - xmin) > 0.00001, "Error: vector x should not have all equal values");
        dbg_assert(degree >= 1, "degree must be >= 1");

        const VectorXf x_norm = ((2*x).array() - (xmax + xmin)) / (xmax - xmin);

        MatrixXf T = MatrixXf::Zero(m, n);
        T.col(0) = VectorXf::Ones(m);
        if (n >= 1) T.col(1) = x_norm;

        for (int j = 2; j < n; ++j) {
            T.col(j) = (2*x_norm).array() * T.col(j-1).array() - T.col(j-2).array();
        }

        // ColPivHouseholderQR<MatrixXf> T_Qr = T.colPivHouseholderQr();
        // dbg_assert(T_Qr.rank() == degree, "");

        HouseholderQR<MatrixXf> T_Qr = T.householderQr();
        dbg_assert(T.colPivHouseholderQr().rank() == degree, "T.colPivHouseholderQr().rank() == degree");

        return chebpoly(T_Qr.solve(y), xmin, xmax);
    }

    VectorXf chebeval(const VectorXf& x, const chebpoly& b, const int degree) {
        dbg_assert(degree >= 1, "degree must be a positive integer");

        const int n = degree;
        const int m = x.rows();
        const float xmax = b.xmax;
        const float xmin = b.xmin;

        dbg_assert(std::abs(xmax - xmin) > 0.00001, "Error: vector x should not have all equal values");
        dbg_assert(degree >= 1, "degree must be >= 1");

        const VectorXf x_norm = ((2*x).array() - (xmax + xmin)) / (xmax - xmin);

        VectorXf y = VectorXf::Zero(m);

        MatrixXf T = MatrixXf::Zero(m, n);
        T.col(0) = VectorXf::Ones(m);
        y += b.coeffs(0) * T.col(0);

        if (n >= 1) {
            T.col(1) = x_norm;
            y += b.coeffs(1) * T.col(1);
        }

        for (int j = 2; j < n; ++j) {
            T.col(j) = (2*x_norm).array() * T.col(j-1).array() - T.col(j-2).array();
            y += b.coeffs(j) * T.col(j);
        }

        return y;
    }
}
