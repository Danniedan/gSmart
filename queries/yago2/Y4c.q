  
SELECT ?name2 WHERE
{
    ?p1 <http://yago-knowledge.org/resource/hasPreferredName> "Cheapo" .
    ?p2 <http://yago-knowledge.org/resource/hasPreferredName> ?name2 .
    ?p1 <http://yago-knowledge.org/resource/isMarriedTo> ?p2 .
    ?p1 <http://yago-knowledge.org/resource/wasBornIn> ?city .
    ?p2 <http://yago-knowledge.org/resource/wasBornIn> ?city .
}
