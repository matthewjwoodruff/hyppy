import qualified Data.Text as Text

data Row = Row { values :: [Text.Text], lineno :: Int } deriving (Show)
data Solution = Solution { row :: Row, objectives :: [Double] } deriving (Show)

-- add a line number to each line
numberedlinesof :: [String] -> Int -> [(String, Int)]
numberedlinesof []      _  = []
numberedlinesof (x:xs)  i  = [(x, i+1)] ++ ( numberedlinesof xs (i+1) )

-- split up a line, retaining its number
sep_row :: Char -> (String, Int) -> Row
sep_row sep (line, linenum) = Row (Text.split (== sep) (Text.pack line)) linenum

-- turn contents into rows
rows :: Char -> String -> [Row]
rows sep contents = map (sep_row sep) (numberedlinesof (lines contents) 0)

-- take objectives and convert to doubles
withoobjectives :: [Int] -> Row -> Solution
withoobjectives indices row = Solution row objs where
    objs = [read (Text.unpack ((values row) !! i)) | i <- indices]

main = do
    contents <- getContents
    print $ (map $ withoobjectives [0]) (rows ' ' contents)

-- vim:ts=4:sw=4:expandtab:ai:number:ruler:
