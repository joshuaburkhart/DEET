require 'test/unit'
require_relative '../lib/MicroArrayHashBuilder'

class MicroArrayHashBuilderTest < Test::Unit::TestCase
    def setup
        #write sample ma_files to disk
        valid_ma_dat1 = <<EOS
gene    M       A       t       pvalue  qvalue  B
F5BTJ3O01A06K0  -5.21883502040817       10.464855458084 -24.3979146617684       3.88889884099112e-41    2.62986784122024e-36    80.8402203948856
CONTIG08318     -4.50800857865597       12.3493124111531        -19.4906231330081       7.85650639514218e-34    2.65648122485745e-29    65.3067831300001
CONTIG08615     -3.60141995355215       10.0230929794592        -18.0331712489824       1.96856524005421e-31    4.43747414528887e-27    60.1100301757823
CONTIG01260     -2.33803021478037       9.93371286896816        -16.473621860716        9.72370858624623e-29    1.64391448286225e-24    54.2290485081661
F5BTJ3O01DFV5K  2.40083486755819        8.45500820725805        20.738433376567 8.48592510614653e-36    5.73860685303159    68.5757063639721
CONTIG15129     2.50319416934623        8.81424858306963        14.9929927590914        4.66804561357376e-26    1.05225528205975    47.8891342577096
F5BTJ3O02J3NEO  2.49613378315458        9.22763092906223        14.4756336362768        4.30904408317478e-25    5.82798212249389    45.7976780181708
F5BTJ3O01B7HBN  2.69803444176811        9.00734381982447        14.202076793316 1.41445842924945e-24    1.5942125212999     44.6766990391772
EOS
        valid_ma_dat2 = <<EOS
gene    M       A       t       pvalue  qvalue  B
CONTIG01260     -2.33803021478037       9.93371286896816        -16.473621860716        9.72370858624623e-29    1.64391448286225    54.2290485081661
CONTIG08318     -4.50800857865597       12.3493124111531        -19.4906231330081       7.85650639514218e-34    2.65648122485745e-29    65.3067831300001
F5BTJ3O01DFV5K  2.40083486755819        8.45500820725805        20.738433376567 8.48592510614653e-36    5.73860685303159e-31    68.5757063639721
F5BTJ3O01A06K0  3.66642347252078       10.464855458084 -17.1404320402277       6.60400966016046e-30    2.23298076634176e-25    56.1630242565358
CONTIG15129     2.50319416934623        8.81424858306963        14.9929927590914        4.66804561357376e-26    1.05225528205975e-21    47.8891342577096
CONTIG08615     -2.95813679393057       10.0230929794592        -14.8120985807975       1.01153535753179e-25    1.71012696382717e-21    47.1621373979199
F5BTJ3O02J3NEO  2.49613378315458        9.22763092906223        14.4756336362768        4.30904408317478e-25    5.82798212249389e-21    45.7976780181708
F5BTJ3O01B7HBN  2.69803444176811        9.00734381982447        14.202076793316 1.41445842924945e-24    1.5942125212999e-20     44.6766990391772
EOS
        invalid_ma_dat = <<EOS
 15 F5BTJ3O01DFV5K  2.40083486755819        8.45500820725805        20.738433376567 8.48592510614653e-36    5.7386    0685303159e-31    68.5757063639721
  16 F5BTJ3O01A06K0  -3.66642347252078       10.464855458084 -17.1404320402277       6.60400966016046e-30    2.2329    8076634176e-25    56.1630242565358
   17 CONTIG15129     2.50319416934623
EOS
        vh1 = File.open("valid_ma_dat1.txt","w")
        vh1.puts(valid_ma_dat1)
        vh1.close
        vh2 = File.open("valid_ma_dat2.txt","w")
        vh2.puts(valid_ma_dat2)
        vh2.close
        ivh = File.open("invalid_ma_dat.txt","w")
        ivh.puts(invalid_ma_dat)
        ivh.close
        @log_filename = "#{Time.now.to_i}.testlog"
        @loghandl = File.open(@log_filename,"w")
    end
    def teardown
        #delete sample ma_files from disk
        %x(rm valid_ma_dat1.txt)
        %x(rm valid_ma_dat2.txt)
        %x(rm invalid_ma_dat.txt)
        @loghandl.close
        %x(rm #{@log_filename})
    end
    def testMakeHash
        #single valid ma_filename
        seq_hash = MicroArrayHashBuilder.makeHash("valid_ma_dat1.txt",@loghandl)
        assert_not_nil(seq_hash)
        assert_equal(Hash,seq_hash.class)
        actual = seq_hash["F5BTJ3O01A06K0"]
        expected = "0"
        assert_equal(expected,actual)
        actual = seq_hash["CONTIG08318"]
        expected = "0"
        assert_equal(expected,actual)
        actual = seq_hash["CONTIG08615"]
        expected = "0"
        assert_equal(expected,actual)
        actual = seq_hash["CONTIG01260"]
        expected = "0"
        assert_equal(expected,actual)
        actual = seq_hash["F5BTJ3O01DFV5K"]
        expected = "X"
        assert_equal(expected,actual)
        actual = seq_hash["CONTIG15129"]
        expected = "X"
        assert_equal(expected,actual)
        actual = seq_hash["F5BTJ3O02J3NEO"]
        expected = "X"
        assert_equal(expected,actual)
        actual = seq_hash["F5BTJ3O01B7HBN"]
        expected = "X"
        assert_equal(expected,actual)
        actual = seq_hash.length
        expected = 8
        assert_equal(expected,actual)

        #multiple valid ma_filenames
        seq_hash = MicroArrayHashBuilder.makeHash("valid_ma_dat1.txt","valid_ma_dat2.txt",@loghandl)
        assert_not_nil(seq_hash)
        assert_equal(Hash,seq_hash.class)
        actual = seq_hash["F5BTJ3O01A06K0"]
        expected = "01"
        assert_equal(expected,actual)
        actual = seq_hash["CONTIG08318"]
        expected = "00"
        assert_equal(expected,actual)
        actual = seq_hash["CONTIG08615"]
        expected = "00"
        assert_equal(expected,actual)
        actual = seq_hash["CONTIG01260"]
        expected = "0X"
        assert_equal(expected,actual)
        actual = seq_hash["F5BTJ3O01DFV5K"]
        expected = "X1"
        assert_equal(expected,actual)
        actual = seq_hash["CONTIG15129"]
        expected = "X1"
        assert_equal(expected,actual)
        actual = seq_hash["F5BTJ3O02J3NEO"]
        expected = "X1"
        assert_equal(expected,actual)
        actual = seq_hash["F5BTJ3O01B7HBN"]
        expected = "X1"
        assert_equal(expected,actual)
        actual = seq_hash.length
        expected = 8
        assert_equal(expected,actual)

        #single invalid ma_filename
        assert_raise ArgumentError do
            MicroArrayHashBuilder.makeHash("an invalid filename",@loghandl)
        end

        #multiple invalid ma_filenames
        assert_raise ArgumentError do
            MicroArrayHashBuilder.makeHash("an invalid filename","another invalid filename",@loghandl)
        end

        #single file with invalid format
        assert_raise ArgumentError do
            MicroArrayHashBuilder.makeHash("invalid_ma_dat.txt",@loghandl)
        end

       #mix valid and invalid file formats
       assert_raise ArgumentError do
           MicroArrayHashBuilder.makeHash("valid_ma_dat1.txt","valid_ma_dat2.txt","invalid_ma_dat.txt",@loghandl)
       end
    end
end
