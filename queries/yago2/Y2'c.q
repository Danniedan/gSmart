SELECT ?GivenName ?FamilyName WHERE{
    ?p <http://yago-knowledge.org/resource/hasGivenName> ?GivenName . 
    ?p <http://yago-knowledge.org/resource/hasFamilyName> ?FamilyName . 
    ?p <http://yago-knowledge.org/resource/wasBornIn> ?city . 
    ?p <http://yago-knowledge.org/resource/hasAcademicAdvisor> <http://yago-knowledge.org/resource/Alfred_Kleiner> .
    <http://yago-knowledge.org/resource/Alfred_Kleiner> <http://yago-knowledge.org/resource/hasGivenName> ?GivenName . 
    ?p <http://yago-knowledge.org/resource/isMarriedTo> ?p2 .
    ?p2 <http://yago-knowledge.org/resource/hasGivenName> ?GivenName . 
}
