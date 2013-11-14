import qualified Data.Text as Text

data Row = Row { values :: [Text.Text]
               , lineno :: Int
               }

data Solution = Solution { row :: Row
                         , objectives :: [Double]
                         }



-- add a line number to each line
numberedlinesof :: [String] -> Int -> [(String, Int)]
numberedlinesof []      _  = []
numberedlinesof (x:xs)  i  = [(x, i+1)] ++ ( numberedlinesof xs (i+1) )

-- split up a line, retaining its number
sep_row :: Char -> (String, Int) -> ([Text.Text], Int)
sep_row sep (line, lineno) = (Text.split (== sep) (Text.pack line), lineno)

-- take a line and convert it to a row
rows :: Char -> String -> [([Text.Text], Int)]
rows sep contents = map (sep_row sep) (numberedlinesof (lines contents) 0)

-- take objectives and convert to doubles
toobjectives :: [Int] -> ([Text.Text], Int) -> ([Text.Text], Int, [Double])
toobjectives indices (row, lineno) = (row, lineno, [read (Text.unpack (row !! i)) | i <- indices])

main = do
    contents <- getContents
    print $ (map $ toobjectives [0]) (rows ' ' contents)

-- vim:ts=4:sw=4:expandtab:ai:number:ruler:
